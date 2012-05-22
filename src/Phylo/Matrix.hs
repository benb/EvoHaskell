module Phylo.Matrix where
import Numeric.LinearAlgebra.LAPACK
import Numeric.LinearAlgebra.Algorithms
import Numeric.LinearAlgebra
import Data.Packed.Matrix
import Data.Packed.Vector
import Data.Packed.ST
import Control.Monad.ST
import Control.Monad
import Debug.Trace 


combineQ :: [Matrix Double] -> Matrix Double
combineQ submats = fromBlocks allSubmats where
                        nonSubmats = zip submats [0..] 
                        zeroMat=buildMatrix rowCount colCount (\(i,j) -> 0)
                        rowCount=rows $ head submats
                        colCount=cols $ head submats
                        tmpSubmats = map (\(mat,i) -> (replicate i zeroMat) ++ (mat: (replicate (len - i -1) zeroMat))) nonSubmats
                        allSubmats = tmpSubmats
                        len = length submats
instance Element Int
interLMat :: Int -> Int -> [[Int]]
interLMat size alphabet = toLists $ fromBlocks allSubmats where
                          diagMat = buildMatrix alphabet alphabet (\(i,j) -> if (i==j) then 1 else 0)
                          zeroMat = buildMatrix alphabet alphabet (\(i,j) -> 0)
                          allSubmats = map (\i -> (replicate i diagMat) ++ (zeroMat : (replicate (size - 1 - i) diagMat))) [0..(size-1)]

intraLMat :: Int -> Int -> [[Int]]
intraLMat size alphabet = toLists $ fromBlocks allSubmats where
                          invDiagMat = buildMatrix alphabet alphabet (\(i,j) -> if (i==j) then 0 else 1)
                          zeroMat = buildMatrix alphabet alphabet (\(i,j) -> 0)
                          allSubmats = map (\i -> (replicate i zeroMat) ++ (invDiagMat : (replicate (size - 1 - i) zeroMat))) [0..(size-1)]

makeQ matS pi = runSTMatrix $ do 
                let rawQ = matS <>(diag pi)
                let rowTots = map (\(row,ind) -> (sumElements row) - (rawQ @@> (ind,ind))) $ zip (toRows (rawQ)) [0..]
                q <- thawMatrix rawQ
                let normRow (i,tot) = writeMatrix q i i (-tot)
                mapM_ (normRow) $ zip [0..(rows rawQ)-1] rowTots
                return q

getRate :: Matrix Double -> Vector Double -> Double
getRate matQ pi = -(sumElements $ zipVectorWith (*) (takeDiag matQ) pi)

setRate rate matQ pi = matMult matQ (rate/oldRate) where
                         oldRate = getRate matQ pi
 
matMult q s = fromRows $ map (mapVector (*s)) $ toRows q

normQ = setRate 1.0

pT :: EigenS -> Double -> Matrix Double
pT (right,lambda,left) t = right <> (diag lambdaT) <> left where
                               lambdaT = mapVector (\x -> (exp (x*t))) lambda

type EigenS = (Matrix Double, Vector Double, Matrix Double)

eigQ :: Matrix Double -> Vector Double -> (Matrix Double, Vector Double, Matrix Double) 
eigQ matQ pi = (u,lambda,u') where
               (myMatA,piRt,piRt') = matA matQ pi
               (lambda,r) = eigS myMatA
               r' = trans r
               u' = r' <> piRt
               u = piRt' <> r

accEigQ matQ pi = distance where
                  (u,lambda,u') = eigQ matQ pi
                  q2' = flatten $ (u <> (diag lambda)) <> u'
                  q' = flatten matQ
                  distance = foldVector (+) 0.0 $ mapVector (\a-> a*a) $  zipVectorWith (-) q2' q'
                  


matA matQ pi = ((piRt <> matQ) <> piRt',piRt,piRt') where
                  piRt = diag (mapVector sqrt pi)
                  invsqrt x = 1/(sqrt x) 
                  piRt' = diag (mapVector invsqrt pi)

symS oldS = runSTMatrix $ do 
            let c = cols oldS
            let r = rows oldS
            s <- thawMatrix oldS
            let ij i j = writeMatrix s i j (oldS @@> (j,i))
            let sym i = mapM_ (ij i) [i..(c-1)]
            mapM_ (sym) [0..(r-1)] 
            return s

fixDiag mat = runSTMatrix $ do 
                let discrepencies = map (foldVector (+) 0) $ toRows mat
                s <- thawMatrix mat
                forM_ (zip discrepencies [0..]) $ \(disc,row) -> do
                        modifyMatrix s row row (\f->f-disc)
                return s
                        

getSwitchingRate mat pi nc = foldVector (+) 0.0 $ zipVectorWith (*) pi (mapVector ( * (-1)) $ takeDiag newmat) where
                                diag = takeDiag newmat
                                newmat = fixDiag $ fromBlocks blocks
                                blocks=(filtBlocks 0) $ toBlocksEvery nc nc mat
                                flat (x:xs) = x ++ (flat xs)
                                flat [] = []
                                empty = (nc><nc) $ repeat 0
                                filtBlocks i (xlist:xs) = (filterList i xlist):(filtBlocks (i+1) xs)
                                filtBlocks i [] = []
                                filterList 0 (x:xs) = empty:xs
                                filterList i (x:xs) = x:(filterList (i-1) xs)
                                filterList i [] = []
                             
