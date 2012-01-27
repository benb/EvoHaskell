module Phylo.Opt (safeGoldenSection,goldenSection,boundedFunction) where
import Debug.Trace
import Data.Maybe
import Data.List
import Numeric.GSL.Minimization

phi = (1+sqrt(5))/2
resphi = 2 - phi

safeGoldenSection tol x1 x3 f = ans where
                                x2 = x1 + (resphi * (x3-x1))
                                f1 = f x1
                                f2 = f x2
                                f3 = f x3
                                ans = if (abs(f3-f2)<tol*10)
                                      then
                                          safeGoldenSection tol x1 (x3/2) f
                                      else 
                                          goldenSection' tol x1 x2 x3 f1 f2 f3 f
                                          

goldenSection :: Double -> Double -> Double -> (Double -> Double) -> Double
goldenSection tol x1 x3 f = goldenSection' tol x1 x2 x3 (f x1) (f x2) (f x3) f where
                                        x2 = x1 + (resphi * (x3-x1))

goldenSection' :: Double -> Double -> Double -> Double -> Double -> Double -> Double -> (Double -> Double) -> Double
goldenSection' tol x1 x2 x3 f1 f2 f3 f = goldenSection'' tol x1 x2 x3 x4 f1 f2 f3 (f x4) f where
                                                x4 = x2 + resphi * (x3 - x2) 

--goldenSection'' tol x1 x2 x3 x4 f1 f2 f3 f4 f | trace ((show [x1,x2,x3,x4]) ++ (show [f1,f2,f3,f4])) False =undefined
goldenSection'' tol x1 x2 x3 x4 f1 f2 f3 f4 f | (abs (x3 - x1) < tol * (abs(x2)+abs(x4))) =  (x3 + x2)/2
                                              | f4 < f2 = goldenSection' tol x2 x4 x3 f2 f4 f3 f
                                              | otherwise = goldenSection' tol x4 x2 x1 f4 f2 f1 f

                                            

boundedFunction :: Double -> [Maybe Double] -> [Maybe Double] -> ([Double] -> Double) -> [Double] -> Double
boundedFunction bad lower upper f vals = if ok
                                         then ans
                                         else bad where
                                         lowerProbs = findIndex (==True) $ map (\(bound,i) -> (fromJust bound) > i)$ filter boundset $ zip lower vals
                                         upperProbs = findIndex (==True) $ map (\(bound,i) -> (fromJust bound) < i)$ filter boundset $ zip upper vals
                                         lowerok = Nothing == lowerProbs
                                         upperok = Nothing == upperProbs
                                         boundset (bound,i) = isJust bound
                                         ok = lowerok && upperok
                                         ans = f vals

