module Data.Chemistry.QC.Psi4
( Psi4Input(..)
, psi4InputParser
, psi4EnergyParser
, psi4GradientParser
, writePsiInput
, calcPsiEnergy
, calcPsiGradient
) where
import System.IO
import System.IO.Unsafe
import System.Process
import Control.Applicative
import Data.Attoparsec.Text.Lazy
import qualified Data.Text as T
import qualified Data.Text.IO as T
import Text.Printf
import qualified Numeric.LinearAlgebra as BLAS
import qualified Data.Chemistry.XYZ as XYZ
import Data.Either.Unwrap


{- ########## -}
{- Data Types -}
{- ########## -}

data Psi4Input = Psi4Input { memory    :: (Int, String)
                           , charge    :: Int
                           , mult      :: Int
                           , molecule  :: XYZ.XYZ
                           , inputopts :: String
                           } deriving Show


{- ###### -}
{- Parser -}
{- ###### -}

psi4InputParser ::  Parser Psi4Input
psi4InputParser = do
    -- parse memory in first line
    skipSpace
    _ <- string $ T.pack "memory"
    skipSpace
    memory_number <- decimal
    skipSpace
    memory_unit <- Data.Attoparsec.Text.Lazy.take 2
    skipSpace
    -- go on to molecule
    _ <- string $ T.pack "molecule"
    skipSpace
    xyz_comment <- manyTill anyChar space
    skipSpace
    _ <- char '{'
    skipSpace
    charge_p <- decimal
    skipSpace
    mult_p <- decimal
    skipSpace
    -- parse molecular coordinates
    mol_content_p <- many' coordLineParser
    _ <- char '}'
    -- parse the set directive for the input options
    inputopts_p <- takeText
    return $ Psi4Input { memory    = (memory_number, T.unpack memory_unit)
                       , charge    = charge_p
                       , mult      = mult_p
                       , molecule  = XYZ.XYZ { XYZ.nAtoms     = length mol_content_p
                                             , XYZ.comment    = xyz_comment
                                             , XYZ.xyzcontent = mol_content_p
                                             }
                       , inputopts = T.unpack inputopts_p
                       }
    where
        coordLineParser :: Parser (String, Double, Double, Double)
        coordLineParser = do
            skipSpace
            element_p <- manyTill anyChar space
            skipSpace
            xp <- double
            skipSpace
            yp <- double
            skipSpace
            zp <- double
            skipSpace
            return $ (element_p, xp, yp, zp)

psi4EnergyParser :: Parser Double
psi4EnergyParser = do
    energies <- many1 psiInterEnergyParser
    return $ last energies
    
    where
        psiInterEnergyParser :: Parser Double
        psiInterEnergyParser = do
            manyTill anyChar (energyString)
            energy <- double
            return $ energy
        
        energyString :: Parser ()
        energyString = do
            _ <- string $ T.pack "Total "
            _ <- manyTill anyChar ((string $ T.pack "Energy") <|> (string $ T.pack "energy"))
            skipSpace
            char '='
            skipSpace
            return ()
        

psi4GradientParser :: Parser (BLAS.Vector Double)
psi4GradientParser = do
    _ <- manyTill anyChar ((string $ T.pack "-Total Gradient:") <|> (string $ T.pack "-Total gradient:"))
    skipSpace
    _ <- string (T.pack "Atom            X                  Y                   Z")
    skipSpace
    _ <- string (T.pack "------   -----------------  -----------------  -----------------")
    skipSpace
    atomgrads <- many' gradientLineParser
    return $ BLAS.fromList $ concat atomgrads
    where
        gradientLineParser :: Parser [Double]
        gradientLineParser = do
            skipSpace
            _ <- (decimal :: Parser Int)
            skipSpace
            gx <- double
            skipSpace
            gy <- double
            skipSpace
            gz <- double
            endOfLine
            return $ [gx, gy, gz]


{- ############ -}
{- IO Functions -}
{- ############ -}

writePsiInput :: Handle -> Psi4Input -> IO ()
writePsiInput handle psi4 = do
    hPrintf handle "%s "      $ "memory"
    hPrintf handle "%d "      $ fst . memory $ psi4
    hPrintf handle "%s\n\n"   $ snd . memory $ psi4
    hPrintf handle "%s "      $ "molecule"
    hPrintf handle "%s "      $ XYZ.comment . molecule $ psi4
    hPrintf handle "%s\n"     $ "{"
    hPrintf handle " %d"      $ charge psi4
    hPrintf handle " %d\n"    $ mult psi4
    mapM_ (\(e, x, y, z) -> do
        hPrintf handle "   %-4s"     e
        hPrintf handle "%+16.8f    " x
        hPrintf handle "%+16.8f    " y
        hPrintf handle "%+16.8f\n"   z) (XYZ.xyzcontent . molecule $ psi4)
    hPrintf handle "%s\n"     $ "}"
    hPutStrLn handle          $ inputopts psi4

calcPsiEnergy :: Psi4Input -> String -> FilePath -> Int -> Double
calcPsiEnergy psiMolecule filePrefix psiPath threads = unsafePerformIO $ do
    let psi4_inputName = (filePrefix ++ "_energy.psi")
        psi4_outputName = (filePrefix ++ "_energy.out")
    
    -- write an input for QC software
    psi4InputHandle <- openFile psi4_inputName WriteMode
    writePsiInput psi4InputHandle psiMolecule
    hClose psi4InputHandle
    
    -- call psi4 and write output
    psiLogging <- readProcess psiPath ["-n", show threads, "-i", psi4_inputName, "-o", psi4_outputName] ""
    
    -- read the output and parse it
    psi4Output <- T.readFile psi4_outputName
    let energy_parsed = parseOnly psi4EnergyParser psi4Output
    
    -- if no parse possible, raise an error beacuse this means the calculation did not work
    energy <- if (isLeft energy_parsed)
       then error $ "calculation did not converge. Please check the psi4 output\n" ++ psiLogging
       else return $ fromRight energy_parsed
    
    -- return the energy from Psi4
    return energy

calcPsiGradient :: Psi4Input -> String -> FilePath -> Int -> BLAS.Vector Double
calcPsiGradient psiMolecule filePrefix psiPath threads = unsafePerformIO $ do
    let psi4_inputName = (filePrefix ++ "_gradient.psi")
        psi4_outputName = (filePrefix ++ "_gradient.out")
    
    -- write an input for QC software
    psi4InputHandle <- openFile psi4_inputName WriteMode
    writePsiInput psi4InputHandle psiMolecule
    hClose psi4InputHandle
    
    -- call psi4 and write output
    psiLogging <- readProcess psiPath ["-n", show threads, "-i", psi4_inputName, "-o", psi4_outputName] ""
    
    -- read the output and parse it
    psi4Output <- T.readFile psi4_outputName
    let gradient_parsed = parseOnly psi4GradientParser psi4Output
    
    -- if no parse possible, raise an error beacuse this means the calculation did not work
    gradient <- if (isLeft gradient_parsed)
       then error $ "calculation did not converge. Please check the psi4 output\n" ++ psiLogging
       else return $ fromRight gradient_parsed
    
    -- return the gradient from Psi4
    return gradient
