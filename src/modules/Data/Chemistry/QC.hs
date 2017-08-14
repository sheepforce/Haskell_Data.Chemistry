{-
implementing the ridge algorithm as described in 
    I. V. Ionova, E. A. Carter, The Journal of Chemical Physics 1993, 98, 6377–6386.
    http://dx.doi.org/10.1063/1.465100

exact implementation is not described well in this paper, so some optimization
procedures are based on my judgement and opinion
-}
module Data.Chemistry.QC
( RelaxMethod(..)
, MicroMethod(..)
, QCSoftware(..)
, Verbosity(..)
, interpolateGeom
, vectorProjection
, steepestDescent
, bfgs
, bfgs2
, makeSaddleNeighbours
, calcDownhillDirection
, relaxPoint
, checkStepJudgement
, list2tuple3
, makePesPoint
, ridge_optimization
) where
import System.IO
import System.IO.Unsafe
import System.Process
import Data.List
import Data.Maybe
import Text.Printf
--import Data.Either
import Data.Either.Unwrap
import Data.Attoparsec.Text.Lazy
import qualified Data.Text as T
import qualified Data.Text.IO as T
import qualified Numeric.LinearAlgebra as BLAS
import qualified Numeric.GSL as GSL
import qualified Data.Chemistry.XYZ as XYZ
import qualified Data.Chemistry.QC.Psi4 as Psi4
import Data.Chemistry.QC.Potential


data QcInMolecule = QcInMolecule { name :: String
                                 , headeropts :: String
                                 , memory :: (Int, String)
                                 , atoms :: [String]
                                 , geom_cart :: BLAS.Vector Double
                                 , multiplicity :: Int
                                 , charge :: Int
                                 , inputopts :: String
                                 }


-- representing a point on the PES with local informations (energy, gradient, hessian)
data PESPoint = PESPoint { loc_geom :: BLAS.Vector Double               -- geometry of this point
                         , loc_energy :: Maybe Double                   -- electronic energy
                         , loc_gradient :: Maybe (BLAS.Vector Double)   -- nuclear gradient
                         , loc_hessian :: Maybe (BLAS.Matrix Double)    -- nuclear hessian at this point
                         } deriving Show

-- representing the three points of one ridge iterations as points on the PES
data Ridge = Ridge { ridge_xs :: PESPoint                               -- the maximum on the connection line
                   , ridge_x0 :: PESPoint                               -- neighbouring point closer to the educt
                   , ridge_x1 :: PESPoint                               -- neighbouring point closer to the product
                   }
{-                                 
methods to use for the relaxation of points x_0' and x_1'
    SaddleInfo = use the gradient and possibly other local informations 
        at the saddle point to relax x_0' and x_1' (equation 3)
    NeighbourInfo = use the gradient and possibly other local informations
        at the the points  x_0' and x_1' directly to relax them (equation 4)
-}
data RelaxMethod = SaddleInfo | NeighbourInfo deriving Eq

-- with which optimizer the algorithm should try to find the maximum on the search path
data MicroMethod = SteepestDescent | BFGS | BFGS2 deriving (Eq, Show)

-- the software used for quantum calculations. Test for a test with the Müller-Brown potential
data QCSoftware = Psi4 | Test deriving Eq

-- how much output should be printed
data Verbosity = Low | Medium | High | Debug deriving Eq

{-
for the relaxation steps of points  x_0' and x_1', choose a uni matrix
as A^T or a full redundant cartesian matrix (would need to be calculated by 
QC software first
-}
--data HessianType = Unit | RedundantCart

-- take two geometries and interpolate on a linear path between them
-- use to make a guess for a good starting point for the optimization to point x*
interpolateGeom :: (BLAS.Vector Double, BLAS.Vector Double) -> Int -> [BLAS.Vector Double]
interpolateGeom (geom1, geom2) nImgs = [(BLAS.vector [p])  * geom2 + (1 - (BLAS.vector [p])) * geom1 | p <- (BLAS.toList $ BLAS.linspace nImgs (0, 1 :: Double))]

vecLength :: BLAS.Vector Double -> Double
vecLength a = sqrt . BLAS.sumElements $ (a * a)

-- project one vector on an other (a on b)
vectorProjection :: BLAS.Vector Double -> BLAS.Vector Double -> BLAS.Vector Double
vectorProjection b_vec a_vec = (BLAS.vector [a1]) * bUnit_vec
    where
        a1 = a_vec BLAS.<.> bUnit_vec
        bUnit_vec = b_vec / (BLAS.fromList [(sqrt . BLAS.sumElements $ (b_vec * b_vec))])

{-
EXAMPLES:
    1) -- optimize to the maximum on the 2D sin surface from start geometry at [-0.5, -0.5] and end geometry at [-0.93, 3.47]
       let start = BLAS.vector [-0.5, -0.5]
       let end = BLAS.vector [-0.93, 3.47]
       let reactionVector = end - start
       steepestDescent 1.0e-6 1000 0.05 1.0e-5 (negate . vec_sin2D) (negate. (vectorProjection reactionVector) . vec_sin2D') (BLAS.vector [-0.59, 0.41])
       (correct result at [-0.706796744008091, 1.4883327295630742])
    2) -- optimize to the maximum on the path between Müller-Brown surface with start at [0.6211317901475866,  2.823674628779631e-2]
       -- and end at [-0.5591208056179865, 1.4407836667987057]
       let start = BLAS.vector [0.6211317901475866,  2.823674628779631e-2]
       let end = BLAS.vector [-0.5591208056179865, 1.4407836667987057]
       let reaction_vec = end - start
       steepestDescent 1.0e-6 1000 0.02 1.0e-5 (negate . muellerBrown (lAK, ak, bk, ck, x0k, y0k)) (negate . (vectorProjection reaction_vec) . muellerBrown' (lAK, ak, bk, ck, x0k, y0k)) end
       (correct result at [-0.19720416827758844, 1.0076355141142914])
-}
-- perform a steepest descent optimization with gradients as implemented in GSL/hmatrix-gsl
-- to make it an ascending algorithm, negate the energy and gradient and project the gradient 
-- on the desired vector
steepestDescent ::
    Double ->                                       -- desired precision of the solution (gradient test)
    Int ->                                          -- maximum number of iterations
    Double ->                                       -- size of first trial step
    Double ->                                       -- tolerance - see GSL manual, recommended 0.1, sets step size
    (BLAS.Vector Double -> Double) ->               -- function to optimize
    (BLAS.Vector Double -> BLAS.Vector Double) ->   -- gradient of the function
    BLAS.Vector Double ->                           -- starting point
    (BLAS.Vector Double, BLAS.Matrix Double)        -- solution vector and optimization path
steepestDescent gradConv maxIter trialStepSize tolerance energy gradient initGeom = 
    GSL.minimizeVD GSL.SteepestDescent gradConv maxIter trialStepSize tolerance energy gradient initGeom

-- BFGS optimization algorithm, also works for the projected maximum when negating energy and gradient and
-- projecting the gradient
bfgs ::
    Double ->                                       -- desired precision of the solution (gradient test)
    Int ->                                          -- maximum number of iterations
    Double ->                                       -- size of first trial step
    Double ->                                       -- tolerance - see GSL manual, recommended 0.1, sets step size
    (BLAS.Vector Double -> Double) ->               -- function to optimize
    (BLAS.Vector Double -> BLAS.Vector Double) ->   -- gradient of the function
    BLAS.Vector Double ->                           -- starting point
    (BLAS.Vector Double, BLAS.Matrix Double)        -- solution vector and optimization path
bfgs gradConv maxIter trialStepSize tolerance energy gradient initGeom = 
    GSL.minimizeVD GSL.VectorBFGS gradConv maxIter trialStepSize tolerance energy gradient initGeom

-- BFGS optimization algorithm, also works for the projected maximum when negating energy and gradient and
-- projecting the gradient
-- RECOMMENDED ALGORITHM
bfgs2 ::
    Double ->                                       -- desired precision of the solution (gradient test)
    Int ->                                          -- maximum number of iterations
    Double ->                                       -- size of first trial step
    Double ->                                       -- tolerance - see GSL manual, recommended 0.1, sets step size
    (BLAS.Vector Double -> Double) ->               -- function to optimize
    (BLAS.Vector Double -> BLAS.Vector Double) ->   -- gradient of the function
    BLAS.Vector Double ->                           -- starting point
    (BLAS.Vector Double, BLAS.Matrix Double)        -- solution vector and optimization path
bfgs2 gradConv maxIter trialStepSize tolerance energy gradient initGeom = 
    GSL.minimizeVD GSL.VectorBFGS gradConv maxIter trialStepSize tolerance energy gradient initGeom
 
{-
create neighbouring points x_0' and x_1' (equation 1)
    x_0' = x* - p
    x_1' = x* + p
    
        distance = side step size
        reaction_vec = the vector between educt and product on which the local maximum was found
        locMax = the coordinates of the local maximum
        
        distance * reaction_vec = p
-}
makeSaddleNeighbours :: Double -> BLAS.Vector Double -> BLAS.Vector Double -> (BLAS.Vector Double, BLAS.Vector Double)
makeSaddleNeighbours distance reaction_vec locMax = (x_0', x_1')
    where
        x_0' = locMax - p
        x_1' = locMax + p
        p = (BLAS.vector [distance]) * reaction_vec

{-
calculate the vectors p_0 and p_1 in (equation 2, 3, 4) with the hessian matrix and the gradient
    hessian = A
    gradient = g (may it be g(x*) or g(x_i'), doesn't matter
-}
calcDownhillDirection :: BLAS.Matrix Double -> BLAS.Vector Double -> BLAS.Vector Double
calcDownhillDirection hessian gradient = negate $ (BLAS.tr' . BLAS.pinv $ hessian) BLAS.#> gradient

{-
calculate new points according to (equation 2)
    stepSize = alpha
    displacement = p
    originalPoint = x_i''
-}
relaxPoint :: Double -> BLAS.Vector Double -> BLAS.Vector Double -> BLAS.Vector Double
relaxPoint stepSize originalPoint displacement = originalPoint + (BLAS.vector [stepSize ]) * displacement

checkStepJudgement :: BLAS.Matrix Double -> BLAS.Vector Double -> (Double, Bool)
checkStepJudgement hessian gradient
    | prod > 0 = (prod, True)
    | prod <= 0 = (prod, False)
    | otherwise = error "checkStepJudgement: can not calculate step criterion"
    where
        grad_mat = BLAS.fromColumns [gradient]
        prod = BLAS.atIndex ((BLAS.tr grad_mat) BLAS.<> (BLAS.pinv hessian) BLAS.<> grad_mat) (0, 0)

list2tuple3 :: [b] -> [(b, b, b)]
list2tuple3 b = [(b !! (ind + 0), b !! (ind + 1), b !! (ind + 2)) | ind <- [0, 3 .. (length b - 1)]]            

-- calculates the energy of a molecule by calling external QC software
-- geometry comes from the vector, not the template
energy_calculator :: QcInMolecule -> QCSoftware -> String -> FilePath -> Int -> BLAS.Vector Double -> Double
energy_calculator molecule software filePrefix programmPath threads geometry
    | software == Psi4 = Psi4.calcPsiEnergy molecule_psi filePrefix programmPath threads
    | software == Test = muellerBrown (lAK, ak, bk, ck, x0k, y0k) geometry
        where
            -- for Psi4
            molecule_psi = Psi4.Psi4Input { Psi4.memory    = memory molecule
                                          , Psi4.charge    = charge molecule
                                          , Psi4.mult      = multiplicity molecule
                                          , Psi4.molecule  = XYZ.XYZ { XYZ.nAtoms     = length . atoms $ molecule
                                                                     , XYZ.comment    = name molecule
                                                                     , XYZ.xyzcontent = 
                                                                         zipWith (\e (x, y, z) -> (e, x, y, z)) (atoms molecule) (list2tuple3 . BLAS.toList $ geometry)
                                                                     }
                                          , Psi4.inputopts = inputopts molecule
                                          }
            -- for Test / Müller-Brown
            lAK = BLAS.vector [-200, -100, -170, 15]
            ak = BLAS.vector [-1, -1, -6.5, 0.7]
            bk = BLAS.vector [0, 0, 11, 0.6]
            ck = BLAS.vector [-10, -10, -6.5, 0.7]
            x0k = BLAS.vector [1, 0, -0.5, -1]
            y0k = BLAS.vector [0, 0.5, 1.5, 1]


-- calculates the energy of a molecule by calling external QC software
gradient_calculator :: QcInMolecule -> QCSoftware -> String -> FilePath -> Int -> BLAS.Vector Double -> BLAS.Vector Double
gradient_calculator molecule software filePrefix programmPath threads geometry
    | software == Psi4 = Psi4.calcPsiGradient molecule_psi filePrefix programmPath threads
    | software == Test = muellerBrown' (lAK, ak, bk, ck, x0k, y0k) geometry
        where
            molecule_psi = Psi4.Psi4Input { Psi4.memory    = memory molecule
                                          , Psi4.charge    = charge molecule
                                          , Psi4.mult      = multiplicity molecule
                                          , Psi4.molecule  = XYZ.XYZ { XYZ.nAtoms     = length . atoms $ molecule
                                                                     , XYZ.comment    = name molecule
                                                                     , XYZ.xyzcontent = 
                                                                         zipWith (\e (x, y, z) -> (e, x, y, z)) (atoms molecule) (list2tuple3 . BLAS.toList $ geometry)
                                                                     }
                                          , Psi4.inputopts = inputopts molecule
                                          }
            lAK = BLAS.vector [-200, -100, -170, 15]
            ak = BLAS.vector [-1, -1, -6.5, 0.7]
            bk = BLAS.vector [0, 0, 11, 0.6]
            ck = BLAS.vector [-10, -10, -6.5, 0.7]
            x0k = BLAS.vector [1, 0, -0.5, -1]
            y0k = BLAS.vector [0, 0.5, 1.5, 1]

-- convert a vector and a template file of a given quantum chemistry input to a QcInMolecule
-- the geometry comes from the Vector, not from the template file
vecTemplate2Qc :: String -> QCSoftware -> [String] -> Int -> Int -> String -> BLAS.Vector Double -> QcInMolecule
vecTemplate2Qc molName software elements mult chrg template geom
    | software == Psi4 = QcInMolecule { name = molName
                                      , headeropts = ""
                                      , memory = Psi4.memory psi4_mol
                                      , atoms = elements
                                      , geom_cart = geom
                                      , multiplicity = mult
                                      , charge = chrg
                                      , inputopts = Psi4.inputopts psi4_mol
                                      }
    | software == Test = QcInMolecule { name = molName
                                      , headeropts = ""
                                      , memory = (1, "Gb")
                                      , atoms = ["X", "Y"]
                                      , geom_cart = geom
                                      , multiplicity = 1
                                      , charge = 0
                                      , inputopts = ""
                                      }
    where
        psi4_mol = fromRight $ parseOnly Psi4.psi4InputParser (T.pack template)

-- helper function to crate PESPoints
--makePesPoint :: BLAS.Vector Double -> Double -> Maybe BLAS.Vector Double -> Maybe BLAS.Matrix Double -> PESPoint
makePesPoint geom energy gradient hessian = PESPoint { loc_geom = geom
                                                     , loc_energy = energy
                                                     , loc_gradient = gradient
                                                     , loc_hessian = hessian
                                                     }

solutionMatrixPrint :: Verbosity -> Handle -> XYZ.XYZ -> BLAS.Matrix Double -> IO ()
solutionMatrixPrint verbosity handle xyzTemplate solMat = do
    let solList = BLAS.toLists solMat
    if (verbosity == Debug)
       then do
           putStrLn                                     $ "  iteration  |    energy                                         "
           putStrLn                                     $ "-------------|--------------                                      "
           mapM_ (\e -> do
               hPrintf handle "%12d | "                 $ (round :: Double -> Int) (e !! 0)
               hPrintf handle "%11.7F \n"               $ negate (e !! 1)
               XYZ.printXYZ stdout (XYZ.vec2Coord xyzTemplate (BLAS.fromList $ drop 2 e))
                 ) solList
    
       else do
           putStrLn                                     $ "  iteration  |    energy                                         "
           putStrLn                                     $ "-------------|--------------                                     "
           mapM_ (\e -> do
               hPrintf handle "%12d | "                 $ (round :: Double -> Int) (e !! 0)
               hPrintf handle "%11.7F \n"               $ negate (e !! 1)
                 ) solList

{-
implementation of the ridge algorithm

the first call supplies start geometries and 
-}
ridge_optimization ::
    -- input values, need to give them at first call
    [String] ->                                     -- elements in the geometries
    BLAS.Vector Double ->                           -- educt geometry
    BLAS.Vector Double ->                           -- product geometry
    Int ->                                          -- multiplicity
    Int ->                                          -- charge
    Int ->                                          -- number of images between product and educt to search for maximum on this path. 1 would be a simple mid point guess
    Int ->                                          -- maximum number of macro iterations
    Int ->                                          -- maximum number of micro iterations
    Double ->                                       -- side step size ||p|| from the paper for building x0' and x1'
    Double ->                                       -- step size for outer steps (alpha)
    Double ->                                       -- step size for inner steps (trust)
    Double ->                                       -- initial step size micro iterations
    Double ->                                       -- inner convergence (gradient)
    (Double, Double) ->                             -- convergence criterium of the gradient (RMS gradient, max single gradient)
    (QCSoftware, String) ->                         -- Quantum chemistry software to use for gradient and possibly hessian calculations and template input file as a single string
    MicroMethod ->                                  -- which optimizer to use for the micro iterations for finding the maximum on the path
    Int ->                                          -- number of threads for the QC software to execute
    Verbosity ->                                    -- how much informations are to be printed?
    -- loop values, given from one iteration to the next to build a history
    Int ->                                          -- macro iteration number !!MUST BE ZERO IN FIRST CALL!!
    [Ridge] ->                                      -- list of maxima and neighbouring points for each macro iteration !!MUST BE EMPTY LIST IN FIRST CALL!!
    [(BLAS.Vector Double, BLAS.Matrix Double)] ->   -- micro iteration info about !!MUST BE EMPTY LIST IN FIRST CALL!!
    -- output
    IO ()
ridge_optimization
    -- input values, given at first iteration and don't change
    elements                                        -- atoms in the correct order in both geometries
    educt                                           -- geometry of the educt
    prod                                            -- geometry of the product
    multiplicity                                    -- multiplicity of the molecule
    charge                                          -- charge of the molecule
    nSearchImages                                   -- for finding the maximum on the path, initially do a calculation of geometries along this path
                                                    -- give a 1 if you want to use the midpoint as a guess and hopefully save 
    maxIterOuter                                    -- maximum number of outer (complete ridge) iterations
    maxIterInner                                    -- maximum number if micro iterations for finding the maximum on the path
    outerSideStepSize                               -- side step size ||p|| from the paper for building x0' and x1'
    outerStepSize                                   -- the alpha from the paper, how large the relaxation shall be
    innerTrust                                      -- trust value of the GSL BFGS or SteepestDescent algorithm
    innerConv                                       -- convergence value of the BFGS
    innerInitStep                                   -- intial search step size for the BFGS or SteepestDescent
    (gradConv_rms, gradConv_max)                    -- convergence of the gradient in Hartree/Bohr
    (software, inputTemplate)                       -- QC software to use and a template file for it (if using the test case, inputTemplate can be anything)
    micromethod                                     -- which optimizer should be used for finding the maximum on the path
    threads
    verbosity                                       -- how much output to print
    -- loop values, cycling through the optimization and updated every step
    nMacroIter                                      -- current macro iteration of the current loop
    ridges                                          -- ridge history (three points: xs, x0, x1) of the optimization
    micro_info                                      -- history of the micro iterations used for finding the maxima
    = do        
        if (nMacroIter == 0)
           then do
               -- check if input is acceptable from a logical point of view and raise error if not"
               if (nSearchImages < 1)
                  then error "number of interpolatet images for searchin the projected maximum is too small"
                  else return ()
               if (maxIterOuter < 1)
                  then error "number of outer iterations is too small"
                  else return ()
               if (maxIterInner < 1)
                  then error "number of inner iterations is too small"
                  else return ()
               if ((gradConv_rms < 0 || gradConv_max < 0) || gradConv_max < gradConv_rms)
                  then error "gradient convergence is non sensical"
                  else return ()
                  
            -- set the microIteration optimization method
            {-
            microOptFunc <- case micromethod of
                                 SteepestDescent -> return steepestDescent
                                 BFGS            -> return bfgs
                                 BFGS2           -> return bfgs2
            -}
                 
               
               -- if high print level is requested the educt and product geometries should be printed
               let educt_temp_xyz = XYZ.XYZ { XYZ.nAtoms = round . (/3) . fromIntegral . length . BLAS.toList $ educt
                                            , XYZ.comment = "educt"
                                            , XYZ.xyzcontent = 
                                                zipWith (\e (x, y, z) -> (e, x, y, z)) elements (list2tuple3 . BLAS.toList $ educt)
                                            }
                   prod_temp_xyz = XYZ.XYZ { XYZ.nAtoms = round . (/3) . fromIntegral . length . BLAS.toList $ prod
                                           , XYZ.comment = "product"
                                           , XYZ.xyzcontent = zipWith (\e (x, y, z) -> (e, x, y, z)) elements (list2tuple3 . BLAS.toList $ prod)
                                           }
                   educt_xyz = XYZ.vec2Coord educt_temp_xyz educt
                   prod_xyz = XYZ.vec2Coord prod_temp_xyz prod
               
               
               putStrLn                                 $ "       ____                                                   "
               putStrLn                                 $ "      < Hi >                                                  "
               putStrLn                                 $ "       ----                                                   "
               putStrLn                                 $ "       \\                                                     "
               putStrLn                                 $ "        \\                                                    "
               putStrLn                                 $ "            __                                               "
               putStrLn                                 $ "           UooU\\.'@@@@@@`.                                  "
               putStrLn                                 $ "           \\__/(@@@@@@@@@@)                                  "
               putStrLn                                 $ "                 (@@@@@@@@)                                   "
               putStrLn                                 $ "                 `YY~~~~YY'                                   "
               putStrLn                                 $ "                  ||    ||                                    "
               putStrLn                                 $ "                                                              "    
               putStrLn                                 $ "          RIDGE OPTIMIZATION                                  "
               putStrLn                                 $ "                                                              "
               putStrLn                                 $ "                                                              "
               if ((verbosity == High || verbosity == Debug )&& software /= Test)
                  then do
                      putStrLn                          $ "your starting geometries are                                  "
                      putStrLn                          $ "    educt:                                                    "
                      putStrLn                          $ "--------------------------------------------------------------"
                      XYZ.printXYZ stdout educt_xyz
                      putStrLn                          $ "--------------------------------------------------------------"
                      putStrLn                          $ "                                                              "
                      putStrLn                          $ "                                                              "
                      putStrLn                          $ "    product:                                                  "
                      putStrLn                          $ "--------------------------------------------------------------"
                      XYZ.printXYZ stdout prod_xyz
                      putStrLn                          $ "--------------------------------------------------------------"
                  else return ()
               putStrLn                                 $ "                                                              "
               putStrLn                                 $ "                                                              "
               putStrLn                                 $ "performing ridge optimization with parameters                 "
               putStrLn                                 $ "    microImages         : " ++ (show nSearchImages)
               putStrLn                                 $ "    max. outer iter     : " ++ (show maxIterOuter)
               putStrLn                                 $ "    max. inner iter     : " ++ (show maxIterInner)
               putStrLn                                 $ "    rms gradient thresh : " ++ (show gradConv_rms)
               putStrLn                                 $ "    max gradient thresh : " ++ (show gradConv_max)
               putStrLn                                 $ "    x* micro algorithm  : " ++ (show micromethod)
               putStrLn                                 $ "                                                              "
               putStrLn                                 $ "                                                              "
               putStrLn                                 $ "                                                              "
               putStrLn                                 $ "                  STARTING RIDGE OPTIMIZATION                 "
               putStrLn                                 $ "                                                              "
               putStrLn                                 $ "MacroIteration : " ++ (show nMacroIter)
               
               -- template XYZ for further printing
               let genericXYZ = XYZ.XYZ { XYZ.nAtoms = length elements
                                        , XYZ.comment = ""
                                        , XYZ.xyzcontent = zipWith (\e (x, y, z) -> (e, x, y, z)) elements $ replicate (length elements) (0.0, 0.0, 0.0)
                                        }
               let programmPath = "/usr/bin/psi4"
               
               -- generate the interpolated geometries and convert them to QC input
               let educt_product_vector = prod - educt
                   maxPath_images = interpolateGeom (educt, prod) nSearchImages
                   pathRidge_QcMol = vecTemplate2Qc ("MacroIter" ++ show nMacroIter ++ "_MaximumPathGRAD") software elements multiplicity charge inputTemplate (head maxPath_images)
                   -- calculate the energies of the images
                   maxPath_energies = map (energy_calculator pathRidge_QcMol software ("MacroIter" ++ show nMacroIter ++ "_MaximumPathGRAD") programmPath threads) maxPath_images
                   -- make this a list of PESPoint
                   maxPath_PES = [makePesPoint (maxPath_images !! ind) (Just (maxPath_energies !! ind)) Nothing Nothing | ind <- [0 .. (length maxPath_images - 1)]]
               
               -- print informations about the interpolated images
               if (verbosity == Debug)
                  then do
                      putStrLn                          $ "    interpolated images                                       "
                      let maxPath_images_xyz = map (XYZ.vec2Coord genericXYZ) maxPath_images
                      mapM_ (XYZ.printXYZ stdout) maxPath_images_xyz
                      putStrLn                          $ "                                                              "
                  else return ()
               if (verbosity == Medium || verbosity == High || verbosity == Debug)
                  then do
                      putStrLn                          $ "    energies of the interpolated images                       "
                      mapM_ (\e -> do
                          putStr "        "
                          print e) $ (map fromJust) . (map loc_energy) $ maxPath_PES
                      putStrLn                          $ "                                                              "
                  else return ()
               
               -- start searching the maximum on the educt -> product vector
               -- find the highest energy of the images
               let maxPath_maxEnergy = maximum maxPath_energies
                   maxPath_maxGeom = maxPath_images !! (fromJust $ findIndex (== maxPath_maxEnergy) maxPath_energies)
               if (verbosity == Medium || verbosity == High || verbosity == Debug)
                  then do
                      putStr                            $ "    maximal energy of the linear path : "
                      print                             $ maxPath_maxEnergy
                      if (verbosity == High || verbosity == Debug)
                         then do
                             let maxPath_maxGeom_xyz = XYZ.vec2Coord genericXYZ maxPath_maxGeom
                             XYZ.printXYZ stdout maxPath_maxGeom_xyz
                         else return ()
                      putStrLn                          $ "                                                                "
                  else return ()
               
               -- start using micro iterations with chosen optimizer to find maximum on the path
               putStrLn                                 $ "    start searching for the projected maximum x*                "
               let projected_optimization = bfgs2 
                                                innerConv
                                                maxIterInner
                                                innerInitStep
                                                innerTrust
                                                (negate . energy_calculator pathRidge_QcMol software ("MacroIter" ++ show nMacroIter ++ "_MaximumPathMicroIter") programmPath threads)
                                                (negate . (vectorProjection educt_product_vector) . gradient_calculator pathRidge_QcMol software ("MacroIter" ++ show nMacroIter ++ "_MaximumPathMicroIter") programmPath threads)
                                                maxPath_maxGeom
                   projected_maximum_geom = fst projected_optimization
                   projected_optimization_history = snd projected_optimization
                   
                   -- calculate the gradient at x*
                   projected_maximum_energy = energy_calculator pathRidge_QcMol software ("MacroIter" ++ show nMacroIter ++ "_xs") programmPath threads projected_maximum_geom
                   projected_maximum_grad = gradient_calculator pathRidge_QcMol software ("MacroIter" ++ show nMacroIter ++ "_xs") programmPath threads projected_maximum_geom
               
               if (verbosity == Medium || verbosity == High || verbosity == Debug)
                  then do
                      if (verbosity == High || verbosity == Debug)
                         then do
                             putStrLn                   $ "        micro iteration optimization for the projected maximum  "
                             solutionMatrixPrint verbosity stdout genericXYZ projected_optimization_history
                         else return ()
                      
                      putStrLn                          $ "                                                                "
                      putStrLn                          $ "        found projected maximum x* with energy                  "
                      putStr                            $ "             "
                      putStrLn                          $ show projected_maximum_energy
                      putStrLn                          $ "                                                                "
                      putStrLn                          $ "        at a geometry of                                        "
                      XYZ.printXYZ stdout (XYZ.vec2Coord genericXYZ projected_maximum_geom)
                      putStrLn                          $ "                                                                "
                      putStrLn                          $ "        and gradient at x*                                      "
                      mapM_ (\(x, y, z) -> do
                          hPrintf stdout "%+16.8f    " x
                          hPrintf stdout "%+16.8f    " y
                          hPrintf stdout "%+16.8f\n"   z) $ list2tuple3 . BLAS.toList $ projected_maximum_grad
                      putStrLn                          $ "                                                                "
                  else return ()
               
               -- make the neighbouring points x0' and x1'
               let x0' = projected_maximum_geom + (BLAS.vector [outerSideStepSize]) * educt_product_vector
                   x1' = projected_maximum_geom - (BLAS.vector [outerSideStepSize]) * educt_product_vector
               
               if (verbosity == High || verbosity == Debug)
                  then do
                      putStrLn                          $ "    generating initially x0' and x1'                            "
                      putStrLn                          $ "        x0' now at                                              "
                      XYZ.printXYZ stdout (XYZ.vec2Coord genericXYZ x0')
                      putStrLn                          $ "                                                                "
                      putStrLn                          $ "        x1' now at                                              "
                      XYZ.printXYZ stdout (XYZ.vec2Coord genericXYZ x1')
                      putStrLn                          $ "                                                                "
                  else return ()
                      
               
               -- calculate the lambda parameter of the location of x*
               let projected_maximum_lambda = (vecLength (projected_maximum_geom - educt)) / (vecLength educt_product_vector)
               if (verbosity == Medium || verbosity == High || verbosity == Debug)
                  then do
                      putStrLn                          $ "                                                                "
                      putStrLn                          $ "        lambda'(x*) is the position of the projected maximum    "
                      putStrLn                          $ "        on the educt (0) -> product (1) path                    "
                      putStrLn                          $ "            lambda'(x*) : " ++ show projected_maximum_lambda
                      putStrLn                          $ "                                                                "
                      putStrLn                          $ "                                                                "
                  else return ()
               
               

               
               
               -- make a ridge
               let ridge_zero = Ridge { ridge_xs = makePesPoint projected_maximum_geom (Just projected_maximum_energy) (Just projected_maximum_grad) Nothing
                                      , ridge_x0 = makePesPoint x0' Nothing Nothing Nothing
                                      , ridge_x1 = makePesPoint x1' Nothing Nothing Nothing
                                      }
               
               putStrLn                                 $ "@    iteraration |  energy(x*)  | rms_grad(x*) | max_grad(x*) | lambda(x*) "
               putStrLn                                 $ "@   -------------|--------------|--------------|--------------|------------"
               printf "@    %11d | %11.7f | %12.7f | %12.7f | %10.6f \n" 
                                                                       nMacroIter
                                                                       projected_maximum_energy
                                                                       (vecLength projected_maximum_grad)
                                                                       (maximum . BLAS.toList $ projected_maximum_grad)
                                                                       projected_maximum_lambda
               
               -- call function again to iterate over ridge applications
               ridge_optimization
                   elements
                   educt
                   prod
                   multiplicity
                   charge
                   nSearchImages
                   maxIterOuter
                   maxIterInner
                   outerSideStepSize
                   outerStepSize
                   innerTrust
                   innerConv
                   innerInitStep
                   (gradConv_rms, gradConv_max)
                   (software, inputTemplate)
                   micromethod
                   threads
                   verbosity
                   --
                   1
                   [ridge_zero]
                   [projected_optimization]
           else do
               -- generate necessary informations first
               let genericXYZ = XYZ.XYZ { XYZ.nAtoms = length elements
                                        , XYZ.comment = ""
                                        , XYZ.xyzcontent = zipWith (\e (x, y, z) -> (e, x, y, z)) elements $ replicate (length elements) (0.0, 0.0, 0.0)
                                        }
                   programmPath = "/usr/bin/psi4"
                   pathRidge_QcMol = vecTemplate2Qc ("MacroIter" ++ show nMacroIter ++ "_MaximumPathGRAD") software elements multiplicity charge inputTemplate educt
               
               
               putStrLn                                 $ "MacroIteration : " ++ (show nMacroIter)
               
               -- relax x0' and x1' (from previous iteration) to x0'' and x1'' (this iteration)
               -- get old xs, x0' and x1' from the previous iteration
               let xs_vec = loc_geom . ridge_xs . head $ ridges
                   x0'_vec = loc_geom . ridge_x0 . head $ ridges
                   x1'_vec = loc_geom . ridge_x1 . head $ ridges
                   x0'_x1'_vec = x1'_vec - x0'_vec
                   xs_lambda = (vecLength (xs_vec - x0'_vec)) / (vecLength x0'_x1'_vec)
               
               if (verbosity == Medium || verbosity == High || verbosity == Debug)
                  then do
                      putStrLn                          $ "    lambda(x*) : " ++ show xs_lambda
                      putStrLn                          $ "                                                                "
                  else return ()
               
               -- if the lambda value of the old ridge is within the simple criteria (0.3 < xs_lambda < 0.7)
               -- then relax by using the gradient form xs to relax x0' and x1'
               -- if not, than calculate the gradients at x0' and x1' and relax it by the individual gradients
               (x0''_vec, x1''_vec) <- if (xs_lambda > 0.3 && xs_lambda < 0.7)
                                          then do
                                              let a_mat = BLAS.ident (length . BLAS.toList $ xs_vec)
                                              let relax_vec = negate $ (BLAS.tr a_mat) BLAS.#> (fromJust . loc_gradient . ridge_xs . head $ ridges) * (BLAS.vector [outerStepSize])
                                                  x0''_vec = x0'_vec + relax_vec
                                                  x1''_vec = x1'_vec + relax_vec
                                              
                                              -- here should be the step judgement by equation (5)
                                              
                                              if (verbosity == Medium || verbosity == High || verbosity == Debug)
                                                 then do
                                                     putStrLn $ "    relaxing x0' and x1' based on the gradient at xs (prev. iter)"
                                                     putStrLn $ "    gradient damping factor is : " ++ (show outerStepSize)
                                                     putStrLn $ "    and the resulting step length is : " ++ (show $ vecLength relax_vec)
                                                     putStrLn $ "                                                                "
                                                 else return ()
                                              
                                              -- new neighbouring points
                                              return (x0''_vec, x1''_vec)
                                          else do
                                              let a_x0'_mat = BLAS.ident (length . BLAS.toList $ xs_vec)
                                                  a_x1'_mat = BLAS.ident (length . BLAS.toList $ xs_vec)
                                              let x0'_gradient = gradient_calculator pathRidge_QcMol software ("MacroIter" ++ show nMacroIter ++ "_x0\'Relax") programmPath threads x0'_vec
                                                  x1'_gradient = gradient_calculator pathRidge_QcMol software ("MacroIter" ++ show nMacroIter ++ "_x0\'Relax") programmPath threads x1'_vec
                                                  x0'_relax_vec = negate $ (BLAS.tr a_x0'_mat) BLAS.#> x0'_gradient * (BLAS.vector [outerStepSize])
                                                  x1'_relax_vec = negate $ (BLAS.tr a_x1'_mat) BLAS.#> x1'_gradient * (BLAS.vector [outerStepSize])
                                                  x0''_vec = x0'_vec + x0'_relax_vec
                                                  x1''_vec = x1'_vec + x1'_relax_vec
                                                  
                                              -- here should be the step judgement by equation (5)
                                              
                                              if (verbosity == Medium || verbosity == High || verbosity == Debug)
                                                 then do
                                                     putStrLn $ "    relaxing x0' and x1' based on the gradients at x0' and x1'"
                                                     putStrLn $ "    gradient damping factor is : " ++ (show outerStepSize)
                                                     putStrLn $ "    and the resulting x0' step length is : " ++ (show $ vecLength x0'_relax_vec)
                                                     putStrLn $ "    and the resulting x1' step length is : " ++ (show $ vecLength x1'_relax_vec)
                                                     putStrLn $ "                                                                "
                                                 else return ()
                                              
                                              -- new neighbouring points
                                              return (x0''_vec, x1''_vec)
               
               -- from x0'' and x1'' search the new x*'
               if (verbosity == Medium || verbosity == High || verbosity == Debug)
                  then do
                      putStrLn                          $ "    x0'' now at                                                 "
                      XYZ.printXYZ stdout (XYZ.vec2Coord genericXYZ x0''_vec)
                      putStrLn                          $ "                                                                "
                      putStrLn                          $ "    x1'' now at                                                 "
                      XYZ.printXYZ stdout (XYZ.vec2Coord genericXYZ x1''_vec)
                  else return ()
               
               -- new interpolation path between them
               let x0''_x1''_interpolImages = interpolateGeom (x0''_vec, x1''_vec) nSearchImages
                   x0''_x1''_imageEnergies = map (energy_calculator pathRidge_QcMol software ("MacroIter" ++ show nMacroIter ++ "_MaximumPathGRAD") programmPath threads) x0''_x1''_interpolImages
                   x0''_x1''_vec = x1''_vec - x0''_vec
               
               -- print informations about the interpolated images
               if (verbosity == Debug)
                  then do
                      putStrLn                          $ "    interpolated images between x0'' and x1''                 "
                      let x0''_x1''_interpolImages_xyz = map (XYZ.vec2Coord genericXYZ) x0''_x1''_interpolImages
                      mapM_ (XYZ.printXYZ stdout) x0''_x1''_interpolImages_xyz
                      putStrLn                          $ "                                                              "
                  else return ()
               if (verbosity == Medium || verbosity == High || verbosity == Debug)
                  then do
                      putStrLn                          $ "    energies of the interpolated images                       "
                      mapM_ (\e -> do
                          putStr "        "
                          print e) $ x0''_x1''_imageEnergies
                      putStrLn                          $ "                                                              "
                  else return ()
               
               -- start searching the maximum on the educt -> product vector
               -- find the highest energy of the images
               let x0''_x1''_maxEnergy = maximum x0''_x1''_imageEnergies
                   x0''_x1''_maxGeom = x0''_x1''_interpolImages !! (fromJust $ findIndex (== x0''_x1''_maxEnergy) x0''_x1''_imageEnergies)
               if (verbosity == Medium || verbosity == High || verbosity == Debug)
                  then do
                      putStr                            $ "    maximal energy of the linear path (x0'' -> x1'') : "
                      print                             $ x0''_x1''_maxEnergy
                      if (verbosity == High || verbosity == Debug)
                         then do
                             let x0''_x1''_maxGeom_xyz = XYZ.vec2Coord genericXYZ x0''_x1''_maxGeom
                             XYZ.printXYZ stdout x0''_x1''_maxGeom_xyz
                         else return ()
                      putStrLn                          $ "                                                                "
                  else return ()
               
               -- start searching the maximum on the new path
               putStrLn                                 $ "    start searching for the projected maximum x*                "
               let projected_optimization = bfgs2 
                                                innerConv
                                                maxIterInner
                                                innerInitStep
                                                innerTrust
                                                (negate . energy_calculator pathRidge_QcMol software ("MacroIter" ++ show nMacroIter ++ "_MaximumPathMicroIter") programmPath threads)
                                                (negate . (vectorProjection x0''_x1''_vec) . gradient_calculator pathRidge_QcMol software ("MacroIter" ++ show nMacroIter ++ "_MaximumPathMicroIter") programmPath threads)
                                                x0''_x1''_maxGeom
                   projected_maximum_geom = fst projected_optimization
                   projected_optimization_history = snd projected_optimization

                   -- calculate the gradient at x*
                   projected_maximum_energy = energy_calculator pathRidge_QcMol software ("MacroIter" ++ show nMacroIter ++ "_xs") programmPath threads projected_maximum_geom
                   projected_maximum_grad = gradient_calculator pathRidge_QcMol software ("MacroIter" ++ show nMacroIter ++ "_xs") programmPath threads projected_maximum_geom
               
               if (verbosity == Medium || verbosity == High || verbosity == Debug)
                  then do
                      if (verbosity == High || verbosity == Debug)
                         then do
                             putStrLn                   $ "        micro iteration optimization for the projected maximum  "
                             solutionMatrixPrint verbosity stdout genericXYZ projected_optimization_history
                         else return ()
                      
                      putStrLn                          $ "                                                                "
                      putStrLn                          $ "        found projected maximum x* with energy                  "
                      putStr                            $ "             "
                      putStrLn                          $ show projected_maximum_energy
                      putStrLn                          $ "                                                                "
                      putStrLn                          $ "        at a geometry of                                        "
                      XYZ.printXYZ stdout (XYZ.vec2Coord genericXYZ projected_maximum_geom)
                      putStrLn                          $ "                                                                "
                      putStrLn                          $ "        and gradient at x*                                      "
                      mapM_ (\(x, y, z) -> do
                          hPrintf stdout "%+16.8f    " x
                          hPrintf stdout "%+16.8f    " y
                          hPrintf stdout "%+16.8f\n"   z) $ list2tuple3 . BLAS.toList $ projected_maximum_grad
                      putStrLn                          $ "                                                                "
                  else return ()
               
               let ridge_new = Ridge { ridge_xs = makePesPoint projected_maximum_geom (Just projected_maximum_energy) (Just projected_maximum_grad) Nothing
                                     , ridge_x0 = makePesPoint x0''_vec Nothing Nothing Nothing
                                     , ridge_x1 = makePesPoint x1''_vec Nothing Nothing Nothing
                                     }
                   
                   xs'_rms_grad = vecLength projected_maximum_grad
                   xs'_max_grad = maximum . BLAS.toList $ projected_maximum_grad
                   
               
               putStrLn                                 $ "@    iteraration |  energy(x*)  | rms_grad(x*) | max_grad(x*) | lambda(x*) "
               putStrLn                                 $ "@   -------------|--------------|--------------|--------------|------------"
               printf "@    %11d | %11.7f | %12.7f | %12.7f | %10.6f \n" 
                                                                       nMacroIter
                                                                       projected_maximum_energy
                                                                       xs'_rms_grad
                                                                       xs'_max_grad
                                                                       xs_lambda
               
               if (xs'_rms_grad < gradConv_rms && xs'_max_grad < gradConv_max)
                  then do
                      putStrLn                          $ "    the final geometry is"
                      XYZ.printXYZ stdout (XYZ.vec2Coord genericXYZ projected_maximum_geom)
                      putStrLn                          $ "                                                                "
                      putStrLn                          $ "       ___________________________________                      "
                      putStrLn                          $ "      < Määäää! Optimization has converged >                    "
                      putStrLn                          $ "       -----------------------------------                      "
                      putStrLn                          $ "       \\                                                       "
                      putStrLn                          $ "        \\                                                      "
                      putStrLn                          $ "            __                                                 "
                      putStrLn                          $ "           UooU\\.'@@@@@@`.                                    "
                      putStrLn                          $ "           \\__/(@@@@@@@@@@)                                    "
                      putStrLn                          $ "                 (@@@@@@@@)                                     "
                      putStrLn                          $ "                 `YY~~~~YY'                                     "
                      putStrLn                          $ "                  ||    ||                                      "
                  else do
                      ridge_optimization
                          elements
                          educt
                          prod
                          multiplicity
                          charge
                          nSearchImages
                          maxIterOuter
                          maxIterInner
                          outerSideStepSize
                          outerStepSize
                          innerTrust
                          innerConv
                          innerInitStep
                          (gradConv_rms, gradConv_max)
                          (software, inputTemplate)
                          micromethod
                          threads
                          verbosity
                          --
                          (nMacroIter + 1)
                          (ridge_new : ridges)
                          (projected_optimization : micro_info)
                      
