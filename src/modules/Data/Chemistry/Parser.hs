module Data.Chemistry.Parser
( xyzParser
, xyzTrajParser
, moldenParser
, nwBasisParser
, gmsBasisParser
, numericalMatrixParser
) where
import           Control.Applicative
import           Data.Attoparsec.Text.Lazy
import           Data.Chemistry.Types
import           Data.Chemistry.Wavefunction
import           Data.Maybe
import qualified Data.Text                   as T
import qualified Numeric.LinearAlgebra       as BLAS
import Numeric.LinearAlgebra ((><))
import Data.List


-- | Make a parser optional, return Nothing if there is no match
maybeOption :: Parser a -> Parser (Maybe a)
maybeOption p = option Nothing (Just <$> p)

-- | give a element and a list. Checks if every element in this list ist the given element
allequal :: (Eq a) => a -> [a] -> Bool
allequal e elist = foldl (\acc x -> if x /= e then False else acc) True elist

--------------------------------------------------------------------------------
-- chemical data formats
--------------------------------------------------------------------------------
{- === -}
{- XYZ -}
{- === -}
xyzParser :: Parser XYZ
xyzParser = do
  skipSpace
  nAtoms_parse <- decimal
  _ <- manyTill anyChar endOfLine
  comment_parse <- manyTill anyChar endOfLine
  --coordinates <- many' xyzCoordLineParser
  coordinates <- count nAtoms_parse xyzCoordLineParser
  -- _ <- endOfLine <|> endOfInput
  return XYZ
    { _xyz_nAtoms = nAtoms_parse
    , _xyz_comment = comment_parse
    , _xyz_xyzcontent = coordinates
    }
    where
      xyzCoordLineParser :: Parser (String,Double,Double,Double)
      xyzCoordLineParser = do
        skipSpace
        element_p <- manyTill anyChar (char ' ')
        skipSpace
        x <- double
        skipSpace
        y <- double
        skipSpace
        z <- double
        skipSpace
        _ <- many' endOfLine
        return (element_p,x,y,z)

xyzTrajParser :: Parser [XYZ]
xyzTrajParser = do
    trajectory <- many' xyzParser
    return trajectory


{- ====== -}
{- Molden -}
{- ====== -}
-- | Parse interesting parts of a Molden file (Atoms, GTO, MOs)
moldenParser :: Parser Molden
moldenParser = do
    -- it all began with "[Molden Format]" ...
    _ <- asciiCI $ T.pack "[Molden Format]"

    -- parse the "[Atoms]" block, containing the geometry informations
    atoms_p <- moldenATOMS

    -- parse the "[GTO]" block
    basfuns_p <- moldenGTO

    -- parse the "[MO]" block
    mmos_p <- moldenMO

    -- return the results
    return Molden
      { _molden_atoms = atoms_p
      , _molden_basfuns = basfuns_p
      , _molden_mos = mmos_p
      }
    where
      moldenTITLE :: Parser String
      moldenTITLE = do
          -- it starts with the "[Title]"
          _ <- manyTill anyChar (asciiCI $ T.pack "[Title]")
          skipSpace
          title_p <- manyTill anyChar endOfLine

          return title_p

      moldenATOMS :: Parser (MoldenUnits, [MoldenCoord])
      moldenATOMS = do
          -- it starts with the "[Atoms]" block
          _ <- manyTill anyChar (asciiCI $ T.pack "[Atoms]")

          -- parse the units. Maybe AU/(AU) or Angs/(Angs)
          skipSpace
          units_p <- do
              units_raw <- manyTill anyChar endOfLine
              if T.isInfixOf (T.toLower . T.pack $ units_raw) (T.pack "(au)")
                 then return Bohr
                 else return Angstrom

          -- the geometry line by line
          skipSpace
          geom <- many1 moldenCoordLineParser

          -- return the results
          return (units_p, geom)

          where
              -- parse a coordinate line in a [Atoms] block
              moldenCoordLineParser :: Parser MoldenCoord
              moldenCoordLineParser = do
                  -- parse the atomic symbol
                  skipSpace
                  element_p <- manyTill anyChar (char ' ')

                  -- parse the number of the atom in the molecule
                  skipSpace
                  natom_p <- decimal

                  -- parse the atomic number (aka number of protons)
                  skipSpace
                  nproton_p <- decimal

                  -- the x coordinate
                  skipSpace
                  x_p <- double

                  -- the y coordinate
                  skipSpace
                  y_p <- double

                  -- the z coordinate
                  skipSpace
                  z_p <- double

                  -- return the results
                  return MoldenCoord
                    { _moldenCoord_element = element_p
                    , _moldenCoord_natom = natom_p
                    , _moldenCoord_nproton = nproton_p
                    , _moldenCoord_coord = (x_p, y_p, z_p)
                    }

      moldenGTO :: Parser [[BasFun]]
      moldenGTO = do
          -- it starts with the "[GTO]" block
          _ <- manyTill anyChar (asciiCI $ T.pack "[GTO]")

          -- possibly there is also a unit here
          _ <- manyTill anyChar endOfLine

          -- parse the basis functions of an atom
          skipSpace
          basfuns_p <- many1 moldenGTOAtomParser

          -- parse the informations about coordinate system
          skipSpace
          _ <- maybeOption (asciiCI $ T.pack "[5D]")
          skipSpace
          _ <- maybeOption (asciiCI $ T.pack "[7F]")
          skipSpace
          _ <- maybeOption (asciiCI $ T.pack "[9G]")
          skipSpace

          -- return the results (basis functions atomwise)
          return basfuns_p
          where
              {- parse the basis functions for a complete atom, for example
                  2  0
                s      3  0
                      13.0100000000         0.0334987264
                       1.9620000000         0.2348008012
                       0.4446000000         0.8136829579
                s      1  0
                       0.1220000000         1.0000000000
                p      1  0
                       0.7270000000         1.0000000000
              -}
              moldenGTOAtomParser :: Parser [BasFun]
              moldenGTOAtomParser = do
                  -- parse the atom number, on which the basis functions are centred
                  skipSpace
                  _ <- (decimal :: Parser Int)

                  -- the zero at the end
                  skipSpace
                  _ <- maybeOption (char '0')
                  skipSpace

                  -- parse the basis functions of the given atom
                  basfunsAtom_p <- many1 moldenGTOAtomBFParser

                  --return the result
                  return basfunsAtom_p


              {- parse a single basis function for a given atom
                s      3  0
                      13.0100000000         0.0334987264
                       1.9620000000         0.2348008012
                       0.4446000000         0.8136829579
              -}
              moldenGTOAtomBFParser :: Parser BasFun
              moldenGTOAtomBFParser = do
                  -- parse the angular momentum of this basis function
                  skipSpace
                  angular_p <- do
                      angular_char <- letter
                      return $ orb2AngMom angular_char

                  -- number of PGTOs to experct
                  skipSpace
                  npgto_p <- decimal

                  -- parse the zero or other number till the endOfLine
                  skipSpace
                  _ <- maybeOption (manyTill anyChar $ char ' ' <|> char '\n')

                  -- parse the PGTOs line by line and make a CGTO of them
                  skipSpace
                  cgto_p <- count npgto_p pgto_and_ContrCoeff

                  return $ BasFun { _basfun_angular = angular_p
                                  , _basfun_radial = cgto_p
                                  }

              -- parse a single pgto and its contraction coefficient
              pgto_and_ContrCoeff :: Parser (PGTO, ContrCoeff)
              pgto_and_ContrCoeff = do
                  -- parse the exponent
                  skipSpace
                  pgto_p <- double

                  -- parse the contraction coefficient
                  skipSpace
                  contrcoeff_p <- double

                  -- return the result
                  return $ (pgto_p, contrcoeff_p)

      moldenMO :: Parser [MoldenMO]
      moldenMO = do
          -- it starts with the "[Atoms]" block
          skipSpace

          _ <- manyTill anyChar (asciiCI $ T.pack "[MO]")

          -- parse the MOs
          skipSpace
          mmos <- many1 singleMMOParser

          -- return the results
          return mmos

          where
              singleMMOParser :: Parser MoldenMO
              singleMMOParser = do
                  -- parse the "Sym" statement
                  skipSpace
                  sym_p <- do
                      sym_temp <- maybeOption symParse
                      if (isNothing sym_temp == True)
                         then return "e"
                         else return $ fromJust sym_temp

                  -- parse the "Ene" statement
                  skipSpace
                  _ <- asciiCI $ T.pack "Ene="
                  skipSpace
                  ene_p <- double

                  -- parse the "Spin" statement
                  skipSpace
                  _ <- asciiCI $ T.pack "Spin="
                  skipSpace
                  spin_p <- do
                      spin_raw <- manyTill anyChar $ (char ' ') <|> (char '\n')
                      if (spin_raw == "Alpha")
                         then return Alpha
                         else return Beta

                  -- parse the "Occup" statement
                  skipSpace
                  _ <- asciiCI $ T.pack "Occup="
                  skipSpace
                  occup_p <- double

                  -- parse the MO coefficients in the AO expansion
                  skipSpace
                  coeffs_p <- do
                      coeffs_raw <- many1 coeff_lineParser
                      return $ BLAS.fromList coeffs_raw

                  -- return the results
                  return MoldenMO
                    { _moldenMO_sym = sym_p
                    , _moldenMO_energy = ene_p
                    , _moldenMO_spin = spin_p
                    , _moldenMO_occup = occup_p
                    , _moldenMO_coeffs = coeffs_p
                    }

              -- wrap the symmetry parser in an optional statement, so it can be optional
              symParse :: Parser String
              symParse = do
                  _ <- asciiCI $ T.pack "Sym="
                  skipSpace
                  sym_p <- manyTill anyChar $ char ' ' <|> char '\n'

                  return $ sym_p


              -- parse a single line of the coefficients in the [MO] block for a single MO
              coeff_lineParser :: Parser Double
              coeff_lineParser = do
                  -- not interested in the first number (label of basis function)
                  skipSpace
                  _ <- (decimal :: Parser Int)

                  -- but parse the mo coefficient
                  skipSpace
                  single_mo_coeff_p <- double

                  -- return the result
                  return single_mo_coeff_p


--------------------------------------------------------------------------------
-- basis sets
--------------------------------------------------------------------------------
-- | NWChem format as obtained from Basis set exchange
-- | parse complete basis set file
nwBasisParser :: Parser [Basis]
nwBasisParser = do
    _ <- manyTill anyChar (string $ T.pack "BASIS \"ao basis\" PRINT")
    --_ <- string $ T.pack "BASIS \"ao basis\" PRINT"
    endOfLine
    basises <- many1 nwSingleBasisParser
    _ <- string $ T.pack "END"

    return basises
    where
      -- parse basis for a single element
      nwSingleBasisParser :: Parser Basis
      nwSingleBasisParser = do
        _ <- string $ T.pack "#BASIS SET: "
        _ <- manyTill anyChar endOfLine
        parsedGauss <- many1 nwConGaussParser
        let elements = map fst parsedGauss
            conGaussians = map snd parsedGauss
        if (allequal (head elements) elements)
          then do
           return Basis
             { element = head elements
             , basFuns = concat conGaussians
             }
          else do
           return Basis
             { element = "X"
             , basFuns = [ConGauss {angMom = 0, expoCoeff_pairs = [(0.0, 0.0)]}]
             }

      {-
      parse o complete block of a contracted gaussian
      SP and other contractions as in the pople basis sets are recognized and flattened
      general contractions are also recognized and flattened
      example:
      C    S
        3047.5249000              0.0018347
         457.3695100              0.0140373

      would give ("C", [ConGauss {angMom = 0, expoCoeff_pairs = [(3047.5249000, 0.0018347), (457.3695100, 0.0140373)]}])

      example:
      C    SP
           7.8682724             -0.1193324              0.0689991

      would give ("C", [ConGauss {angMom = 0, expoCoeff_pairs = [(7.8682724, -0.1193324)]}, ConGauss {angMom = 1, expoCoeff_pairs = [(7.8682724, 0.0689991)]}])

      example:
      C    S
           0.1687144              1.0000000              2.0000000

      would give ("C", [ConGauss {angMom = 1, expoCoeff_pairs = [(0.1687144, 1.0)], ConGauss {angMom = 0, expoCoeff_pairs = [(0.1687144,1.0)]}])
      -}
      nwConGaussParser :: Parser (String, [ConGauss])
      nwConGaussParser = do
        atom <- manyTill anyChar (char ' ')
        _ <- many1 (char ' ')
        angMom_raw <- manyTill anyChar (char ' ' <|> char '\n')
        expoCoeffs_pairs_raw <- many1 nwConGaussLineParser
        conGauss_list <- if (length angMom_raw == 1)
          then do
            if ((length . snd $ head expoCoeffs_pairs_raw) == 1)
              then do
                let expoList = map fst expoCoeffs_pairs_raw
                    coeffList = concat . (map snd) $ expoCoeffs_pairs_raw
                return [ ConGauss
                           { angMom = orb2AngMom $ head angMom_raw
                           , expoCoeff_pairs = zip expoList coeffList
                           }
                       ]
              else do
                let expoList = map fst expoCoeffs_pairs_raw
                    manyCoeffList = map snd expoCoeffs_pairs_raw
                return [ ConGauss
                           { angMom = orb2AngMom $ head angMom_raw
                           , expoCoeff_pairs = zip expoList (map (!! ind) manyCoeffList)
                           }
                       | ind <- [0..(length (head manyCoeffList) - 1)]
                       ]
          else do
           let angMoms = map orb2AngMom angMom_raw
               expoList = map fst expoCoeffs_pairs_raw
               manyCoeffList = map snd expoCoeffs_pairs_raw
           return [ ConGauss
                      { angMom = angMoms !! ind
                      , expoCoeff_pairs = zip expoList (map (!! ind) manyCoeffList)
                      }
                  | ind <- [0..(length angMom_raw -1)]
                  ]
        return (atom, conGauss_list)

      {-
      parse single line of gaussians in NWChem style
      example:
           exponent1             coefficient1a              coefficient1b

      would give (exponent, [coefficient1a, coefficient1b])
      -}
      nwConGaussLineParser :: Parser (Double, [Double])
      nwConGaussLineParser = do

        _ <- many' (char ' ')
        expo <- double
        _ <- many' (char ' ')
        coefflist <- many1 coeffparser
        endOfLine
        let pairs = (expo, coefflist)
        return pairs
        where
          coeffparser :: Parser Double
          coeffparser = do
          coeff <- double
          _ <- many' (char ' ')
          return coeff


-- | parse a basis set for GAMESS-US
gmsBasisParser :: Parser [Basis]
gmsBasisParser = do
  _ <- manyTill anyChar (string $ T.pack "$DATA")
  endOfLine
  basises <- many1 gmsSingleBasisParser
  _ <- string $ T.pack "$END"
  return basises
  where
    gmsSingleBasisParser :: Parser Basis
    gmsSingleBasisParser = do
      skipSpace
      fullElemName <- manyTill anyChar endOfLine
      parsedGauss <- many1 gmsConGaussParser
      return Basis
        { element = fullElemName
        , basFuns = concat parsedGauss
        }

    gmsConGaussParser :: Parser [ConGauss]
    gmsConGaussParser = do
        skipSpace
        angMom_raw <- anyChar
        skipSpace
        _ <- (decimal :: Parser Int)
        _ <- many' (char ' ')
        endOfLine
        expoCoeffs_pairs_raw <- many1 gmsConGaussLineParser
        conGauss_list <- if (angMom_raw /= 'L')
          then do
            if (length . snd $ head expoCoeffs_pairs_raw) == 1
              then do
               let expoList = map fst expoCoeffs_pairs_raw
                   coeffList = concat . (map snd) $ expoCoeffs_pairs_raw
               return [ ConGauss
                          { angMom = orb2AngMom angMom_raw
                          , expoCoeff_pairs = zip expoList coeffList
                          }
                      ]
              else do
               let expoList = map fst expoCoeffs_pairs_raw
                   manyCoeffList = map snd expoCoeffs_pairs_raw
               return [ ConGauss
                          { angMom = orb2AngMom angMom_raw
                          , expoCoeff_pairs = zip expoList (map (!! ind) manyCoeffList)
                          }
                      | ind <- [0..(length (head manyCoeffList) - 1)]
                      ]
          else do
            let angMoms = [0, 1]
                expoList = map fst expoCoeffs_pairs_raw
                manyCoeffList = map snd expoCoeffs_pairs_raw
            return [ ConGauss
                       { angMom = angMoms !! ind
                       , expoCoeff_pairs = zip expoList (map (!! ind) manyCoeffList)
                       }
                   | ind <- [0, 1]
                   ]
        return conGauss_list

    gmsConGaussLineParser :: Parser (Double, [Double])
    gmsConGaussLineParser = do
      skipSpace
      _ <- (decimal :: Parser Int)
      skipSpace
      expo <- double
      skipSpace
      coefflist <- many1 coeffparser
      _ <- many' (char ' ')
      endOfLine
      let pairs = (expo, coefflist)
      return pairs
      where
        coeffparser :: Parser Double
        coeffparser = do
          coeff <- double
          _ <- many' (char ' ')
          return coeff


--------------------------------------------------------------------------------
-- Simple data formats
--------------------------------------------------------------------------------
-- | A simple parser for many fixed column data,
numericalMatrixParser :: [Char] -> Parser (Maybe (BLAS.Matrix Double))
numericalMatrixParser commentChar = do
  _ <- many' commentLineParser
  numLines <- many1 numLineParser
  if (length . nub $ [length (numLines !! i) | i <- [0 .. length numLines - 1]]) /= 1
    then return Nothing
    else return $ Just $ ((length numLines) >< (length . head $ numLines)) (concat numLines)
  where
    commentLineParser = do
      skipSpace
      _ <- satisfy (`elem` commentChar)
      commentLine <- manyTill anyChar endOfLine
      return commentLine

    numLineParser = do
      _ <- many' (char ' ' <|> char '\t')
      numLine <- many1 $ do
        n <- double
        _ <- many1 (satisfy (`elem` ";, \t"))
        return n
      endOfLine
      return numLine
