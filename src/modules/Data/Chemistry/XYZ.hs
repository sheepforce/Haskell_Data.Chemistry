module Data.Chemistry.XYZ
( XYZ(..)
, xyzParser
, xyzTrajParser
, getNAtoms
, getElements
, getCoords
, coord2Mat
, mat2Coord
, interpolate
, align
, findNeighboursOfType_Distance
, distance
, angle
) where
import           Data.Chemistry.Parser
import           Data.Chemistry.Types
import           Lens.Micro.Platform
import qualified Numeric.LinearAlgebra as BLAS


getNAtoms :: XYZ -> Int
getNAtoms a = _xyz_nAtoms a

getElement :: (String,Double,Double,Double) -> String
getElement (e, _, _, _) = e

getElements :: XYZ -> [String]
getElements a = map getElement $ a ^. xyz_xyzcontent

getCoord :: (String,Double,Double,Double) -> (Double,Double,Double)
getCoord (_, x, y, z) = (x, y, z)

getCoords :: XYZ -> [(Double,Double,Double)]
getCoords a = map getCoord $ a ^. xyz_xyzcontent

-- | Generate HMatrix Matrix of coordinates
coord2Mat :: XYZ -> BLAS.Matrix Double
coord2Mat a = (_xyz_nAtoms a) BLAS.>< 3 $ concat . map (\(x, y, z) -> [x, y, z]) $ getCoords a

-- | Fill XYZ with coordinates from a matrix (replacing the old ones)
mat2Coord :: XYZ -> BLAS.Matrix Double -> XYZ
mat2Coord a mat = XYZ
  { _xyz_nAtoms = (a ^. xyz_nAtoms)
  , _xyz_comment = (a ^. xyz_comment)
  , _xyz_xyzcontent = content
  }
    where
        content = zipWith (\e (x, y, z) -> (e, x, y, z)) (getElements a) coords
        coords = (listList2TupleList . BLAS.toLists) mat

        listList2TupleList :: [[Double]] -> [(Double, Double, Double)]
        listList2TupleList listlist = map (\[x,y,z] -> (x,y,z)) listlist


--------------------------------------------------------------------------------
-- Working with XYZ files
--------------------------------------------------------------------------------
-- | interpolate two geometries in cartesian coordinates
interpolate :: Int -> XYZ -> XYZ -> [XYZ]
interpolate nImages mol1 mol2 = map (mat2Coord mol1) interpol_mats
    where
        interpol_mats = [mol1_coordMat * ((1 BLAS.>< 1) [(1.0 - p)]) + mol2_coordMat * ((1 BLAS.>< 1) [p]) | p <- (BLAS.toList $ BLAS.linspace nImages (0, 1 :: Double))]
        mol1_coordMat = coord2Mat mol1
        mol2_coordMat = coord2Mat mol2

-- | align molecule by three atoms
-- | first goes to the origin
-- | second goes to the z axis
-- | third goes to the xz plane
align :: (Int, Int, Int) -> XYZ -> XYZ
align (atom1, atom2, atom3) mol = mat2Coord mol mol_aligned
    where
        -- transform to matrix for easy algebra
        mol_mat = coord2Mat mol

        -- translate so that atom 1 is in the origin
        atom1_coords = mol_mat BLAS.! atom1
        mol_1toOrigin = BLAS.fromRows $ [mol_mat BLAS.! ind - atom1_coords | ind <- [0 .. (BLAS.rows mol_mat - 1)]]

        -- rotate that atom 2 on Z-axis
        -- first bring atom 2 to x = 0 by rotation around z
        -- get coordinates from atom 2 and then the necessary angle for this rotation
        atom2_coords = mol_1toOrigin BLAS.! atom2
        x2 = atom2_coords BLAS.! 0
        y2 = atom2_coords BLAS.! 1
        rotAngleZ_2toX0 = (-2) * atan((y2 - sqrt(x2**2 + y2**2)) / x2)
        mat_rotZ_2toX0 = (3 BLAS.>< 3) [ cos rotAngleZ_2toX0, (-1) * sin rotAngleZ_2toX0, 0
                                       , sin rotAngleZ_2toX0, cos rotAngleZ_2toX0       , 0
                                       , 0                  , 0                         , 1]
        mol_2toX0_vecs = map (mat_rotZ_2toX0 BLAS.#>) $ BLAS.toRows mol_1toOrigin
        mol_2toX0 = BLAS.fromRows mol_2toX0_vecs

        -- second bring atom 2 to y = 0 by rotation around x
        atom2_coords_s = mol_2toX0 BLAS.! atom2
        y2_s = atom2_coords_s BLAS.! 1
        z2_s = atom2_coords_s BLAS.! 2
        rotAngleX_2toY0 = (-2) * atan((z2_s - sqrt(y2_s**2 + z2_s**2)) / y2_s)
        mat_rotX_2toY0 = (3 BLAS.>< 3) [ 1, 0, 0
                                       , 0, cos rotAngleX_2toY0, (-1) * sin rotAngleX_2toY0
                                       , 0, sin rotAngleX_2toY0, cos rotAngleX_2toY0]
        mol_2toY0_vecs = map (mat_rotX_2toY0 BLAS.#>) $ BLAS.toRows mol_2toX0
        mol_2toY0 = BLAS.fromRows mol_2toY0_vecs

        -- now that atom 2 is on the z axis, bring atom3 to XZ plane by rotation around z
        atom3_coords = mol_2toY0 BLAS.! atom3
        x3 = atom3_coords BLAS.! 0
        y3 = atom3_coords BLAS.! 1
        rotAngleZ_3toY0 = (-2) * atan((y3 - sqrt(x3**2 + y3**2)) / x3)
        mat_rotZ_3toY0 = (3 BLAS.>< 3) [ cos rotAngleZ_3toY0 , (-1) * sin rotAngleZ_3toY0 , 0
                                       , sin rotAngleZ_3toY0 , cos rotAngleZ_3toY0        , 0
                                       , 0                   , 0                          , 1]
        mol_3toY0_vecs = map (mat_rotZ_3toY0 BLAS.#>) $ BLAS.toRows mol_2toY0
        mol_3toY0 = BLAS.fromRows mol_3toY0_vecs

        -- give the final rotated structure
        mol_aligned = mol_3toY0

-- | get a XYZ file with all atoms of element "t" within distance "d" of the
-- | "a"th atom in the xyz file "xyz"
findNeighboursOfType_Distance :: Int -> String -> Double -> XYZ -> Maybe XYZ
findNeighboursOfType_Distance a t d xyz
  | a < length content && length allOfTypeWithinDistance > 0 = Just XYZ
    { _xyz_nAtoms = length allOfTypeWithinDistance
    , _xyz_comment = "all " ++ t ++ " within " ++ show d ++ " Angstrom of atom " ++ show a
    , _xyz_xyzcontent = allOfTypeWithinDistance
    }
  | otherwise = Nothing
  where
    content = xyz ^. xyz_xyzcontent
    aCoords = content !! a
    allOfTypeWithinDistance =
      filter (\i -> distance (r3Vec2hmVec . getCoord $ aCoords) (r3Vec2hmVec . getCoord $ i) <= d) .
      filter (\(e, _, _, _) -> e == t) $
      content


{-
started here to implement sorting function for wannier centres. Nevertheless,
this is already in TRAVIS

-- | Sort atoms in trajectory frames according to their distance. The atoms of
-- | frame 2 are reshuffled to have minimal distance to the same type in frame 1
sortTrajAtomsByDistance :: [XYZ] -> BLAS.Matrix Double --  [XYZ]
sortTrajAtomsByDistance t =
  f2fDistMat (coordsOfAllFrames !! 0) (coordsOfAllFrames !! 1)
  where
    nAtoms = getNAtoms . head $ t
    coordsOfAllFrames = map coord2Mat t

    f2fDistMat fA fB =
      BLAS.fromLists
      [
        [ distance (fARows !! a) (fBRows !! b)
        | b <- [ 0 .. nAtoms - 1 ]
        ]
      | a <- [ 0 .. nAtoms - 1 ]
      ]
      where
        fARows = BLAS.toRows fA
        fBRows = BLAS.toRows fB
    shuffleMap =
-}
--------------------------------------------------------------------------------
-- Math helper functions
--------------------------------------------------------------------------------
r3Vec2hmVec :: (Double, Double, Double) -> BLAS.Vector Double
r3Vec2hmVec (a_x, a_y, a_z) = BLAS.fromList [a_x, a_y, a_z]

vecLength :: BLAS.Vector Double -> Double
vecLength = sqrt . sum . map (^2) . BLAS.toList

-- | calculate the distance between two points
distance :: BLAS.Vector Double -> BLAS.Vector Double -> Double
distance a b = vecLength (a - b)

-- | calculate angle between two vectors
angle :: BLAS.Vector Double -> BLAS.Vector Double -> Double
angle a b =
  acos $
  (a BLAS.<.> b) /
  ((vecLength a) * (vecLength b))
