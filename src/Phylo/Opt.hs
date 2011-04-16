module Phylo.Opt (goldenSection) where
import Debug.Trace

phi = (1+sqrt(5))/2
resphi = 2 - phi

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

                                            

