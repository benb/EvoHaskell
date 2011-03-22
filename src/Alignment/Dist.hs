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
hom0Dist = zeroDist numberifyBasic
homDist = genDist numberifyBasic
homGapDist = genDist numberifyGap
homTreeDist t = tupDist (numberifyGapTree t)

isPermutation :: ListAlignment -> ListAlignment -> Bool
isPermutation a b | (names a) /= (names b) = False
isPermutation a b = (dropGaps a) == (dropGaps b)


genDist :: (ListAlignment -> [[(Int)]]) -> ListAlignment -> ListAlignment -> (Int,Int)
genDist f = tupDist (tupify f Just)

zeroDist f = tupDist (tupify f (\x->Nothing))
                        
tupify :: (ListAlignment -> [[(Int)]]) -> (Int->Maybe Int) -> (ListAlignment -> [[(Int,Maybe Int)]])
tupify f gapHandler = fmap (map (map toTup)) f where 
        toTup i = if (i < 0) then
                        (i,gapHandler i)
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
--diffIn' (x:xs) (y:ys) ans | trace ((show x) ++ (show y)) False = undefined
--First, skip gaps on left side of xs and ys
diffIn' (((x1,m1),(x2,xx2)):xs) y (i,j) | x1<0 || m1/=Nothing = i `seq` j `seq` diffIn' xs y (i,j)
diffIn' x (((y1,m1),y2):ys) (i,j) | y1<0 || m1/=Nothing  = i `seq` j `seq` diffIn' x ys (i,j)

--Metric 0, ignore gaps
diffIn' ((x1,(x2i,x2g)):xs) ((y1,(y2i,y2g)):ys) (i,j) | (x2i<0 && y2i>0 && x2g==Nothing) || (x2i>0 && y2i<0 && y2g==Nothing)  = i `seq` j `seq` diffIn' xs ys  (i+1,j+1)
--Same
diffIn' ((x1,x2):xs) ((y1,y2):ys) (i,j) | x2==y2  = i `seq` j `seq` diffIn' xs ys  (i+2,j)
--Different
diffIn' ((x1,x2):xs) ((y1,y2):ys) (i,j)  = i `seq` j `seq` diffIn' xs ys (i+2,j+2) 
--diffIn' xs ys ans | trace ((show xs) ++ (show ys)) False = undefined
diffIn' [] [] t = t
