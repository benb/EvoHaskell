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

makeQ matS pi = runSTMatrix $ do 
                let rawQ = matS <>(diag pi)
                let rowTots = map (sumElements) $ toRows (rawQ)
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
                        

