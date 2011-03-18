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


tupDistSeq :: (Eq a, Show a) => [(Int,Maybe a)] -> [(Int,Maybe a)] -> [[(Int,Maybe a)]] ->  [[(Int,Maybe a)]] -> (Int,Int) -> (Int,Int)
tupDistSeq seqA seqB (seqA2:xs) (seqB2:ys) ans = tupDistSeq seqA seqB xs ys $! (diffIn' (zip seqA seqA2) (zip seqB seqB2) ans)
tupDistSeq seqA seqB [] [] ans = ans

tupDist' ::  (Eq a, Show a) => [[(Int,Maybe a)]] ->  [[(Int,Maybe a)]] -> [[(Int,Maybe a)]] ->  [[(Int,Maybe a)]] -> (Int,Int) -> (Int,Int)
--tupDist' (x:xs) (y:ys) headx heady t | trace (show t) False = undefined
--tupDist' (x:xs) (y:ys) headx heady (0,0) | False  = undefined
tupDist' (x:xs) (y:ys) headx heady (i,j) = i `seq` j `seq` tupDist' xs ys (x:headx) (y:heady) (tupDistSeq x y xs ys (tupDistSeq x y headx heady (i,j)))
tupDist' [] [] headx heady t = t 


tupDist :: (Eq a,Show a) => (ListAlignment -> [[(Int,Maybe a)]]) -> ListAlignment -> ListAlignment -> (Int,Int)
tupDist numF aln1 aln2 =  tupDist' num1 num2 [] [] (0,0) where
                                num1 = numF aln1
                                num2 = numF aln2
        

diffIn' :: (Eq a,Show a)=> [((Int,Maybe a),(Int,Maybe a))] -> [((Int,Maybe a),(Int,Maybe a))] -> (Int,Int) -> (Int,Int)
--Gap
diffIn' (((x1,Just f),(x2,xx2)):xs) y (i,j) = i `seq` j `seq` diffIn' xs y (i,j)
diffIn' x (((y1,Just f),y2):ys) (i,j) = i `seq` j `seq` diffIn' x ys (i,j)
--Same
diffIn' ((x1,x2):xs) ((y1,y2):ys) (i,j) | x2==y2  = i `seq` j `seq` diffIn' xs ys  (i+1,j)
--Different
diffIn' ((x1,x2):xs) ((y1,y2):ys) (i,j)  = i `seq` j `seq` diffIn' xs ys (i+1,j+1) 
diffIn' [] [] t = t
