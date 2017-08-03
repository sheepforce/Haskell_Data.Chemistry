module Data.Chemistry.BasisSet
( Basis(..)
, ConGauss(..)
, nwBasisParser
, printBagelBasisList
) where
import qualified Data.Text as T
import Data.Attoparsec.Text.Lazy
import System.IO
import Control.Applicative
import Text.Printf


{- ################# -}
{- define Data Types -}
{- ################# -}

-- Define data representing a finite gaussian basis
data Basis = Basis { element  :: String
                   , basFuns :: [ConGauss]
                   } deriving Show

-- a single contracted gaussian made of primitve gaussians
-- (Double, Double) contains exponents and coefficients
data ConGauss = ConGauss { angMom          :: Int
                         , expoCoeff_pairs :: [(Double, Double)]
                         } deriving Show


{- ###### -}
{- Parser -}
{- ###### -}

-- NWChem format as obtained from Basis set exchange
-- parse complete basis set file
nwBasisParser :: Parser [Basis]
nwBasisParser = do
    _ <- string $ T.pack "BASIS \"ao basis\" PRINT"
    endOfLine
    basises <- many1 nwSingleBasisParser
    _ <- string $ T.pack "END"
    return $ basises

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
           return $ Basis { element = head elements
                          , basFuns = concat conGaussians
                          }
       else do
           return $ Basis { element = "X"
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
                                     return $ [ ConGauss { angMom = orb2AngMom $ head angMom_raw
                                                         , expoCoeff_pairs = zip expoList coeffList
                                                         }
                                              ]
                                 else do
                                     let expoList = map fst expoCoeffs_pairs_raw
                                         manyCoeffList = map snd expoCoeffs_pairs_raw
                                     return $ [ ConGauss { angMom = orb2AngMom $ head angMom_raw
                                                         , expoCoeff_pairs = zip expoList (map (!! ind) manyCoeffList)
                                                         }
                                                | ind <- [0..(length (head manyCoeffList) - 1)]
                                              ]
                          else do
                              let angMoms = map orb2AngMom angMom_raw
                                  expoList = map fst expoCoeffs_pairs_raw
                                  manyCoeffList = map snd expoCoeffs_pairs_raw
                              return $ [ ConGauss { angMom = angMoms !! ind
                                                  , expoCoeff_pairs = zip expoList (map (!! ind) manyCoeffList)
                                                  }
                                         | ind <- [0..(length angMom_raw -1)]
                                       ]
    return $ (atom, conGauss_list)

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
    return $ pairs
    where
        coeffparser :: Parser Double
        coeffparser = do
            coeff <- double
            _ <- many' (char ' ') 
            return $ coeff


{- ############### -}
{- other functions -}
{- ############### -}
-- convert from orbital type to quantum number
orb2AngMom :: Char -> Int
orb2AngMom orb
    | orb == 's' || orb == 'S' = 0
    | orb == 'p' || orb == 'P' = 1
    | orb == 'd' || orb == 'D' = 2
    | orb == 'f' || orb == 'f' = 3
    | orb == 'g' || orb == 'g' = 4
    | orb == 'h' || orb == 'h' = 5
    | orb == 'i' || orb == 'i' = 6
    | orb == 'j' || orb == 'j' = 7
    | otherwise = error "not supported angular momentum"

-- angular momentum to orbital
angMom2Orb :: Int -> Char
angMom2Orb angmom
    | angmom == 0 = 's'
    | angmom == 1 = 'p'
    | angmom == 2 = 'd'
    | angmom == 3 = 'f'
    | angmom == 4 = 'g'
    | angmom == 5 = 'h'
    | angmom == 6 = 'i'
    | angmom == 7 = 'j'
    | otherwise = error "not supported angular momentum"

-- give a element and a list. Checks if every element in this list ist the given element    
allequal :: (Eq a) => a -> [a] -> Bool
allequal e elist = foldl (\acc x -> if x /= e then False else acc) True elist


{- ############ -}
{- IO functions -}
{- ############ -}
-- give a list of a atomic basises and print them in a json format, used by Bagel
printBagelBasisList :: Handle -> [Basis] -> IO ()
printBagelBasisList handle basises = do
    hPrintf handle "%s\n" $ "{"
    printBasisList {-handle-} basises
    hPrintf handle "\n%s\n" $ "}"
    where
        -- print a list of basises to json, but omit the top level braces
        printBasisList :: {-Handle ->-} [Basis] -> IO ()
        printBasisList {-handle-} [] = return ()
        printBasisList {-handle-} [a] = printBasis {-handle-} a
        printBasisList {-handle-} (a:b) = do
            printBasis {-handle-} a
            hPrintf handle "%s\n" $ ","
            printBasisList {-handle-} b
        
        -- print the basis of a single element but avoid braces of higher levels
        printBasis :: {-Handle ->-} Basis -> IO ()
        printBasis {-handle-} basis = do
            hPrintf handle "  %s\n"  $ "\"" ++ element basis ++ "\" : ["
            printConGaussList {-handle-} $ basFuns basis
            hPrintf handle "\n  %s"    $ "]"
        
        -- print the contracted gaussians but not the element or braces of higher levels
        printConGaussList :: {-Handle ->-} [ConGauss] -> IO ()
        printConGaussList {-handle-} []    = return ()
        printConGaussList {-handle-} [a]   = printConGauss {-handle-} a
        printConGaussList {-handle-} (a:b) = do
            printConGauss {-handle-} a
            hPrintf handle "%s\n" $ ","
            printConGaussList {-handle-} b
        
        -- print a single contracted gaussian in json/Bagel format but not braces of higher levels
        printConGauss :: {-Handle ->-} ConGauss -> IO ()
        printConGauss {-handle-} congauss  = do
            -- print opening brace
            hPrintf handle "    %s\n"     $ "{"
            -- print angular momentum
            hPrintf handle "      %s\n"   $ "\"angular\" : \"" ++ [(angMom2Orb $ angMom congauss)] ++ "\","
            -- print exponents of the primitve gaussians
            hPrintf handle "      %s"     $ "\"prim\" : ["
            printDoubleList {-handle-} (map fst $ expoCoeff_pairs congauss)
            hPrintf handle "%s\n"         $ "],"
            -- print contraction coefficients
            hPrintf handle "      %s"     $ "\"cont\" : [["
            printDoubleList {-handle-} (map snd $ expoCoeff_pairs congauss)
            hPrintf handle "%s\n"         $ "]]"
            -- print closing brace
            hPrintf handle "    %s"       $ "}"
        
        -- print a list of doubles with 7 decimals precision as it would look in haskell and json
        printDoubleList :: {-Handle ->-} [Double] -> IO ()
        printDoubleList {-handle-} []    = return ()
        printDoubleList {-handle-} [a]   = hPrintf handle "%.7F" a
        printDoubleList {-handle-} (a:b) = do
            hPrintf handle "%.7F, " a
            printDoubleList {-handle-} b
        

