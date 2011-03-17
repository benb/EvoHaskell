module Alignment.Dist where
import Alignment
import Debug.Trace

--
--gapEvent :: Node -> ListAlignment -> [[Maybe Split]]
--gapEvent tree (ListAlignment names seqs cols) = gapEvent' tree names cols
--
--gapEvent' :: Node -> [Name] -> [Column] -> [[Maybe Split]]
--gapEvent' (Tree l r) names cols = if (contained leftNames gapNames)  where
--                                             gapNames = map (\x -> fst x) $ filter (\x -> (snd x=='-')) $ zip names cols 
--                                             leftNames = names l
--                                             rightNames = names r
--                                             contained x y = contained' x x y
--                                             contained' full (x:[]) (y:[]) = x==y
--                                             contained' full (x:[]) (y:ys) = x==y || contained full full ys
--                                             contained' full (x:xs) (y:ys) = x==y || contained full xs (y:ys)
--
homDist = genDist numberifyBasic
homGapDist = genDist numberifyGap
homTreeDist t = tupDist (numberifyGapTree t)



genDist :: (ListAlignment -> [[(Int)]]) -> ListAlignment -> ListAlignment -> (Int,Int)
genDist f = tupDist (tupify f)
                        
tupify :: (ListAlignment -> [[(Int)]]) -> (ListAlignment -> [[(Int,Maybe Int)]])
tupify f = fmap (map (map toTup)) f where 
        toTup i = if (i < 0) then
                        (i,Just i)
                  else 
                        (i,Nothing)

tupDist :: (Eq a,Show a) => (ListAlignment -> [[(Int,Maybe a)]]) -> ListAlignment -> ListAlignment -> (Int,Int)
tupDist numF aln1 aln2 = answer (diff3 pairs1 pairs2) where
        
        pairs1 = pairs (numF aln1)
        pairs2 = pairs (numF aln2)

        answer (x,y) = (y,x) -- for some reason I reversed the convention inside tupDist

        pairs :: [[b]] -> [[[(b,b)]]]
        pairs [] = []
        pairs list = pairs' [] list
        pairs' :: [[b]] -> [[b]] -> [[[(b,b)]]]
        pairs' head (x:xs) = (pairsXY (zip (repeat x) (reverse head)) ++ pairsXY (zip (repeat x) xs)) : (pairs' (x:head) xs) 
        pairs' head [] = [] --(pairsXY (zip (repeat x) (reverse head)) ++ pairsXY (zip (repeat x) xs)) : (pairs' x:head xs)

        pairsXY []  = []
        pairsXY ((a,b):xs) = zip a b : pairsXY xs

        addT (a,b) (c,d) = (a+c,b+d)

        diff2 :: (Eq a,Show a)=> [((Int,Maybe a),(Int,Maybe a))] -> [((Int,Maybe a),(Int,Maybe a))] -> (Int,Int)
        --skip gaps
        diff2 (((x1,Just f),(x2,xx2)):xs) y = diff2 xs y
        diff2 x (((y1,Just f),y2):ys) = diff2 x ys
        --Same
--        diff2 ((x1,x2):xs) ((y1,y2):ys) | x1 == y1 && x2==y2 && trace ("Same " ++ (show x1) ++ " " ++ (show y1) ++ " " ++ (show x2) ++ " " ++ (show y2)) False = undefined
        diff2 ((x1,x2):xs) ((y1,y2):ys) | x1 == y1 && x2==y2 = addT (diff2 xs ys) (1,0)
        --Different
--        diff2 ((x1,x2):xs) ((y1,y2):ys) | x1 == y1 && trace ("Diff " ++ (show x1) ++ " " ++ (show y1) ++ " " ++ (show x2) ++ " " ++ (show y2)) False = undefined
        diff2 ((x1,x2):xs) ((y1,y2):ys) | x1 == y1 = addT (diff2 xs ys) (1,1)
        diff2 [] [] = (0,0)

        diff :: (Eq a,Show a) => [[((Int,Maybe a),(Int,Maybe a))]] -> [[((Int,Maybe a),(Int,Maybe a))]] -> (Int,Int)
        diff [] [] = (0,0)
        diff (x:xs) (y:ys) = addT (diff2 x y) (diff xs ys)

        
        diff3 :: (Eq a,Show a) => [[[((Int,Maybe a),(Int,Maybe a))]]] -> [[[((Int,Maybe a),(Int,Maybe a))]]] -> (Int,Int)
        diff3 [] [] = (0,0)
        diff3 (x:xs) (y:ys) = addT (diff x y) (diff3 xs ys)
