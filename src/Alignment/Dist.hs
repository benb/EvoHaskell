module Alignment.Dist where
import Alignment

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
homTreeDist t = tupDist (bodgeTupify $ numberifyGapTree t)



genDist :: (ListAlignment -> [[(Int)]]) -> ListAlignment -> ListAlignment -> Double
genDist f = tupDist (tupify f)
                        
tupify :: (ListAlignment -> [[(Int)]]) -> (ListAlignment -> [[(Int,Maybe Int)]])
tupify f = fmap (map (map toTup)) f where 
        toTup i = if (i < 0) then
                        (i,Just i)
                  else 
                        (i,Nothing)

--Not sure what best API is - TODO fix this!
bodgeTupify :: (ListAlignment -> [[Either Int a]]) -> (ListAlignment -> [[(Int,Maybe a)]])
bodgeTupify f = fmap (map (map toTup)) f where
        toTup (Right a) = (-1,Just a)
        toTup (Left i) = (i,Nothing)



tupDist :: Eq a => (ListAlignment -> [[(Int,Maybe a)]]) -> ListAlignment -> ListAlignment -> Double
tupDist numF aln1 aln2 = answer (diff3 pairs1 pairs2) where
        
        pairs1 = pairs (numF aln1)
        pairs2 = pairs (numF aln2)

        answer :: (Int,Int) -> Double
        answer (numPairs,numDiffs) = (fromIntegral numDiffs) / (fromIntegral numPairs)
        
        pairs :: [[b]] -> [[[(b,b)]]]
        pairs [] = []
        pairs (x:[]) = []
        pairs (x:xs) = pairsXY (zip (repeat x) xs) : pairs (xs)

        pairsXY []  = []
        pairsXY ((a,b):xs) = zip a b : pairsXY xs

        addT (a,b) (c,d) = (a+c,b+d)

        diff2 :: Eq a=> [((Int,Maybe a),(Int,Maybe a))] -> [((Int,Maybe a),(Int,Maybe a))] -> (Int,Int)
        --skip gaps
        diff2 (((x1,Just f),(x2,xx2)):xs) y = diff2 xs y
        diff2 x (((y1,Just f),y2):ys) = diff2 x ys
        --Same
        diff2 ((x1,x2):xs) ((y1,y2):ys) | x1 == y1 && x2==y2 = addT (diff2 xs ys) (1,0)
        --Different
        diff2 ((x1,x2):xs) ((y1,y2):ys) | x1 == y1 = addT (diff2 xs ys) (1,1)
        diff2 [] [] = (0,0)

        diff :: Eq a => [[((Int,Maybe a),(Int,Maybe a))]] -> [[((Int,Maybe a),(Int,Maybe a))]] -> (Int,Int)
        diff [] [] = (0,0)
        diff (x:xs) (y:ys) = addT (diff2 x y) (diff xs ys)

        
        diff3 :: Eq a => [[[((Int,Maybe a),(Int,Maybe a))]]] -> [[[((Int,Maybe a),(Int,Maybe a))]]] -> (Int,Int)
        diff3 [] [] = (0,0)
        diff3 (x:xs) (y:ys) = addT (diff x y) (diff3 xs ys)
