{-# LANGUAGE FlexibleInstances #-} 
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
homDist = labDist numberifyBasic
homGapDist = labDist numberifyGap
homTreeDist t = labDist (numberifyGapTree t)

isPermutation :: ListAlignment -> ListAlignment -> Bool
isPermutation a b | (names a) /= (names b) = False
isPermutation a b = (dropGaps a) == (dropGaps b)


genDist :: (ListAlignment -> [[(Int)]]) -> ListAlignment -> ListAlignment -> (Int,Int)
genDist f = labDist f

zeroDist f = labDistTrig (siteLabel f (\x->Nothing))
                        
siteLabel :: (ListAlignment -> [[(Int)]]) -> (Int->Maybe Int) -> (ListAlignment -> [[(Int,Maybe Int)]])
siteLabel f gapHandler = fmap (map (map toLabel)) f where 
        toLabel i = if (i < 0) then
                        (i,gapHandler i)
                    else 
                        (i,Nothing)

class (Eq a, Show a) => SiteLabel a where 
  isGap :: a -> Bool

instance SiteLabel Int where
  isGap a = a<0


instance (Integral a, Eq b, Show b) => SiteLabel (a,b) where
  isGap (a,b) = a<0


labDistSeq :: (SiteLabel a) => [a] -> [a] -> [[a]] ->  [[a]] -> (Int,Int) -> (Int,Int)
labDistSeq seqA seqB (seqA2:xs) (seqB2:ys) ans = labDistSeq seqA seqB xs ys $! (diffIn' (zip seqA seqA2) (zip seqB seqB2) ans)
labDistSeq seqA seqB [] [] ans = ans

labDist' ::  (SiteLabel a) => [[(a)]] ->  [[(a)]] -> [[(a)]] ->  [[(a)]] -> (Int,Int) -> (Int,Int)
--labDist' (x:xs) (y:ys) headx heady t | trace (show t) False = undefined
--labDist' (x:xs) (y:ys) headx heady (0,0) | False  = undefined
labDist' (x:xs) (y:ys) headx heady (i,j) = i `seq` j `seq` labDist' xs ys (x:headx) (y:heady) (labDistSeq x y xs ys (labDistSeq x y headx heady (i,j)))
labDist' [] [] headx heady t = t 

labDistTrig' ::  (SiteLabel a) => [[(a)]] ->  [[(a)]] -> (Int,Int) -> (Int,Int)
--labDist' (x:xs) (y:ys) headx heady t | trace (show t) False = undefined
--labDist' (x:xs) (y:ys) headx heady (0,0) | False  = undefined
labDistTrig' (x:xs) (y:ys) (i,j) | trace ((show x)++(show y)) False = undefined 
labDistTrig' (x:xs) (y:ys) (i,j) = i `seq` j `seq` labDistTrig' xs ys (labDistSeq x y xs ys (i,j))
labDistTrig' [] [] t = t 


labDistTrig numF aln1 aln2 =  labDistTrig' num1 num2 (0,0) where
                                num1 = numF aln1
                                num2 = numF aln2
        

labDist :: (SiteLabel a) => (ListAlignment -> [[(a)]]) -> ListAlignment -> ListAlignment -> (Int,Int)
labDist numF aln1 aln2 =  labDist' num1 num2 [] [] (0,0) where
                                num1 = numF aln1
                                num2 = numF aln2
        

diffIn' :: (SiteLabel a)=> [(a,a)] -> [(a,a)] -> (Int,Int) -> (Int,Int)
diffIn' (x:xs) (y:ys) ans | trace ((show x) ++ (show y) ++ (show ans)) False = undefined
--First, skip gaps on left side of xs and ys
diffIn' ((x1,x2):xs) y (i,j) | (isGap x1) = i `seq` j `seq` diffIn' xs y (i,j)
diffIn' x ((y1,y2):ys) (i,j) | (isGap y1) = i `seq` j `seq` diffIn' x ys (i,j)

--Metric 0, ignore gaps
--diffIn' ((x1,(x2i,x2g)):xs) ((y1,(y2i,y2g)):ys) (i,j) | (x2i<0 && y2i>0 && x2g==Nothing) || (x2i>0 && y2i<0 && y2g==Nothing)  = i `seq` j `seq` diffIn' xs ys  (i+1,j+1)
--Same
diffIn' ((x1,x2):xs) ((y1,y2):ys) (i,j) | x2==y2  = i `seq` j `seq` diffIn' xs ys  (i+2,j)
--Different
diffIn' ((x1,x2):xs) ((y1,y2):ys) (i,j)  = i `seq` j `seq` diffIn' xs ys (i+2,j+2) 
--diffIn' xs ys ans | trace ((show xs) ++ (show ys)) False = undefined
diffIn' [] [] t = t
