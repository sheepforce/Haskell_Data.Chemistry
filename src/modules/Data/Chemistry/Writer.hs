module Data.Chemistry.Writer
( write_XYZ
, write_Molden
, write_BagelBasis
) where
import           Data.Chemistry.Types
import           Data.Chemistry.Wavefunction
import           Lens.Micro.Platform
import qualified Numeric.LinearAlgebra       as BLAS
import           System.IO
import           Text.Printf

--------------------------------------------------------------------------------
-- Writing Molden files
--------------------------------------------------------------------------------
write_Molden :: Molden -> String
write_Molden molden =
  -- begin with oblgigatory header
  "[Molden Format]" ++ "\n" ++

  -- title is not of intereset and therefore empty
  "[Title]" ++ "\n" ++
  "\n" ++

  -- the [Atoms] section with the coordinates
  "[Atoms] " ++
  printf "%-4s  "
    ( if fst (molden ^. molden_atoms) == Bohr
      then "AU"
      else "Angs"
    ) ++
  "\n" ++
  concat
    ( map (\(n, a) -> printf "%-3s    " (a ^. moldenCoord_element) ++
                      printf "%4d     " n ++
                      printf "%4d     " (a ^. moldenCoord_nproton) ++
                      printf "%+14.8E    " (a ^. moldenCoord_coord . _1) ++
                      printf "%+14.8E    " (a ^. moldenCoord_coord . _2) ++
                      printf "%+14.8E    " (a ^. moldenCoord_coord . _3) ++
                      "\n"
          ) numbAtoms
    ) ++ "\n" ++

  -- the [GTO] section with the basis functions
  "[GTO]" ++ "\n" ++
  concat
    ( map (\(n, bfs) -> printf "  %4d  " n ++ "0\n" ++
                        basFunsWriter bfs ++ "\n"
          ) numbBasFuns
    ) ++
  "[5D]\n" ++
  -- "[7F]\n" ++
  "[9G]\n" ++
  "\n" ++

  -- the [MO] section with the molecular orbitals
  "[MO]\n" ++
  concat (map moWriter orbs)
  where
    mAtoms = snd $ molden ^. molden_atoms
    atomIndRange = [ 1 .. length mAtoms ]
    numbAtoms = zip atomIndRange mAtoms
    numbBasFuns = zip atomIndRange (molden ^. molden_basfuns)
    orbs = molden ^. molden_mos

    -- write an MO
    moWriter :: MoldenMO -> String
    moWriter mo =
      "Sym= " ++ (mo ^. moldenMO_sym) ++ "\n" ++
      "Ene= " ++ printf "%20.10E" (mo ^. moldenMO_energy) ++ "\n" ++
      "Spin= " ++ show (mo ^. moldenMO_spin) ++ "\n" ++
      "Occup= " ++ printf "%8.6F" (mo ^. moldenMO_occup) ++ "\n" ++
      concat
        ( map (\(n, c) -> printf "%5d        " n ++
                          printf "%20.16F    " c ++ "\n"
              ) numbCoeffs
        )
      where
        coeffIndRange = [1 .. ( length . BLAS.toList $ (mo ^. moldenMO_coeffs))]
        numbCoeffs = zip coeffIndRange (BLAS.toList $ mo ^. moldenMO_coeffs)



    -- Write all basis functions of an atom. Results for example in
    {-
    s  1  1.0
            0.9059000000            1.0000000000
    s  1  1.0
            0.1285000000            1.0000000000
    p  5  1.0
           18.7100000000            0.0140309984
            4.1330000000            0.0868659900
            1.2000000000            0.2902159662
            0.3827000000            0.5010079423
            0.1209000000            0.3434059602
    p  1  1.0
            0.3827000000            1.0000000000
    p  1  1.0
            0.1209000000            1.0000000000
    d  1  1.0
            1.0970000000            1.0000000000
    d  1  1.0
            0.3180000000            1.0000000000
    f  1  1.0
            0.7610000000            1.0000000000
    -}
    basFunsWriter :: [BasFun] -> String
    basFunsWriter bfs =
      concat . map basFunWriter $ bfs

    -- Write a single basis function. A basis function contains many nPrimitives
    -- Results for example in
    {-
    s  10 1.0
         8236.0000000000         0.0005309999
         1235.0000000000         0.0041079989
          280.8000000000         0.0210869945
           79.2700000000         0.0818529787
           25.5900000000         0.2348169389
            8.9970000000         0.4344008868
            3.3190000000         0.3461289099
            0.9059000000         0.0393779897
            0.3643000000        -0.0089829977
            0.1285000000         0.0023849994
    -}
    basFunWriter :: BasFun -> String
    basFunWriter bf =
      [angmomString] ++ "  " ++ show nPrimitives ++ "  1.0\n" ++
      concat
        ( map (\(e, c) -> printf "%20.10F    " e ++
                          printf "%20.10F\n"   c
              ) radialPrimitives
        )
      where
        angmom = bf ^. basfun_angular
        angmomString = angMom2Orb angmom
        nPrimitives = length (bf ^. basfun_radial)
        radialPrimitives = bf ^. basfun_radial


--------------------------------------------------------------------------------
-- Writing XYZ files
--------------------------------------------------------------------------------
-- | Print XYZ formated to a handle
-- | Using mapM_ over a list enables printing a trajectory
write_XYZ :: XYZ -> String
write_XYZ xyz =
  show (xyz ^. xyz_nAtoms) ++ "\n" ++
  show (xyz ^. xyz_comment) ++ "\n" ++
  concat
    ( map (\(e, x, y, z) -> printf "%-4s" e ++
                            printf "%+16.8f    " x ++
                            printf "%+16.8f    " y ++
                            printf "%+16.8f\n" z
          ) (xyz ^. xyz_xyzcontent)
    )

--------------------------------------------------------------------------------
-- Bagels JSON basis set format
--------------------------------------------------------------------------------
-- give a list of a atomic basises and print them in a json format, used by Bagel
write_BagelBasis :: Handle -> [Basis] -> IO ()
write_BagelBasis handle basises = do
  hPrintf handle "%s\n" $ "{"
  printBasisList basises
  hPrintf handle "\n%s\n" $ "}"
  where
    -- print a list of basises to json, but omit the top level braces
    printBasisList :: {-Handle ->-} [Basis] -> IO ()
    printBasisList [] = return ()
    printBasisList [a] = printBasis a
    printBasisList (a:b) = do
      printBasis a
      hPrintf handle "%s\n" $ ","
      printBasisList b

    -- print the basis of a single element but avoid braces of higher levels
    printBasis :: {-Handle ->-} Basis -> IO ()
    printBasis basis = do
      hPrintf handle "  %s\n"  $ "\"" ++ element basis ++ "\" : ["
      printConGaussList $ basFuns basis
      hPrintf handle "\n  %s"    $ "]"

    -- print the contracted gaussians but not the element or braces of higher levels
    printConGaussList :: {-Handle ->-} [ConGauss] -> IO ()
    printConGaussList []    = return ()
    printConGaussList [a]   = printConGauss a
    printConGaussList (a:b) = do
      printConGauss a
      hPrintf handle "%s\n" $ ","
      printConGaussList b

    -- print a single contracted gaussian in json/Bagel format but not braces of higher levels
    printConGauss :: {-Handle ->-} ConGauss -> IO ()
    printConGauss congauss  = do
      -- print opening brace
      hPrintf handle "    %s\n"     $ "{"
      -- print angular momentum
      hPrintf handle "      %s\n"   $ "\"angular\" : \"" ++ [(angMom2Orb $ angMom congauss)] ++ "\","
      -- print exponents of the primitve gaussians
      hPrintf handle "      %s"     $ "\"prim\" : ["
      printDoubleList (map fst $ expoCoeff_pairs congauss)
      hPrintf handle "%s\n"         $ "],"
      -- print contraction coefficients
      hPrintf handle "      %s"     $ "\"cont\" : [["
      printDoubleList (map snd $ expoCoeff_pairs congauss)
      hPrintf handle "%s\n"         $ "]]"
      -- print closing brace
      hPrintf handle "    %s"       $ "}"

    -- print a list of doubles with 7 decimals precision as it would look in haskell and json
    printDoubleList :: {-Handle ->-} [Double] -> IO ()
    printDoubleList []    = return ()
    printDoubleList [a]   = hPrintf handle "%.7F" a
    printDoubleList (a:b) = do
      hPrintf handle "%.7F, " a
      printDoubleList b
