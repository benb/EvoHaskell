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

zeroDist f = labDistTrig f
                        
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


instance (Integral a, Eq b, Show b,Ord a) => SiteLabel (a,b) where
  isGap (a,b) = a<0



labDistTrig' ::  (SiteLabel a, Ord a) => [[(a)]] ->  [[(a)]] -> (Int,Int) -> (Int,Int)
---labDistTrig' (x:xs) (y:ys) (i,j) | trace ((show x)++ " , " ++(show y) ++ " " ++ (show (i,j))) False = undefined 
labDistTrig' (x:xs) (y:ys) (i,j) = i `seq` j `seq` labDistTrig' xs ys (labDistSeq diffSSP x y xs ys (i,j))
labDistTrig' [] [] t = t 


labDistTrig numF aln1 aln2 =  labDistTrig' num1 num2 (0,0) where
                                num1 = numF aln1
                                num2 = numF aln2
        

-- |'labDist' computes the distance between two alignments after labelling
-- it takes a labelling function and two alignments
-- and returns a a tuple of (denonimator,numerator), i.e. distance is
-- snd/fst
labDist :: (SiteLabel a) => (ListAlignment -> [[(a)]]) -> ListAlignment -> ListAlignment -> (Int,Int)
labDist numF aln1 aln2 | trace "Fast dist" False  = undefined
labDist numF aln1 aln2 =  labDist' num1 num2 [] [] (0,0) where
                                num1 = numF aln1
                                num2 = numF aln2

-- |takes four lists (x, y, remaining_x and remaining_y) and a distance tuple
-- and increments the distance tuple
labDist' ::  (SiteLabel a) => [[(a)]] ->  [[(a)]] -> [[(a)]] ->  [[(a)]] -> (Int,Int) -> (Int,Int)
--labDist' (x:xs) (y:ys) headx heady t | trace (show t) False = undefined
--labDist' (x:xs) (y:ys) headx heady (0,0) | False  = undefined
labDist' (x:xs) (y:ys) headx heady (i,j) = i `seq` j `seq` labDist' xs ys (x:headx) (y:heady) (labDistSeq diffIn x y xs ys (labDistSeq diffIn x y headx heady (i,j)))
labDist' [] [] headx heady t = t 

-- | increments the tuple (final arg) with the distance between
-- two lists of labels and each of the corresponding lists-of-lists of labels
labDistSeq :: ([a]->[a]->[a]->[a]->(Int,Int)->(Int,Int))-> [a] -> [a] -> [[a]] ->  [[a]] -> (Int,Int) -> (Int,Int)
labDistSeq f seqA seqB (seqA2:xs) (seqB2:ys) ans = labDistSeq f seqA seqB xs ys $! (f seqA seqA2 seqB seqB2 ans)
labDistSeq f seqA seqB [] [] ans = ans
       
-- | compute distance for pairs of labels 
diffIn :: (SiteLabel a) => [a] -> [a] -> [a] -> [a] -> (Int,Int) -> (Int,Int)
diffIn x1s x2s y1s y2s t = addT t $ diffIn' x1s x2s y1s y2s 
--diffIn (x:xs) (y:ys) ans | trace ((show x) ++ (show y) ++ (show ans)) False = undefined
--First, skip gaps on left side of xs and ys
diffIn' (x1:x1s) (x2:x2s) y1s y2s  | (isGap x1) = diffIn' x1s x2s y1s y2s
diffIn' x1s x2s (y1:y1s) (y2:y2s)  | (isGap y1) = diffIn' x1s x2s y1s y2s 
--Same
diffIn' (x1:x1s) (x2:x2s) (y1:y1s) (y2:y2s) | x2==y2  = addT (2,0) $ diffIn' x1s x2s y1s y2s 
--Different
                                            | otherwise = addT (2,2) $ diffIn' x1s x2s y1s y2s 
diffIn' [] [] [] [] = (0,0)

addT (i,j) (i2,j2) = (i+i2,j+j2)



-- | compute distance for pairs of labels for metric 0 (SSP)
diffSSP :: (SiteLabel a, Ord a) => [a] -> [a] -> [a] -> [a] -> (Int,Int) -> (Int,Int)
--diffSSP (x:xs) (y:ys) (i,j) | trace ("diffSSP" ++ (show x) ++ " , " ++  (show y) ++ " " ++ (show (i,j))) False = undefined
diffSSP (a:x1s) (b:x2s) y1s y2s (i,j) | isGap a || isGap b = diffSSP x1s x2s y1s y2s (i,j)
diffSSP x1s x2s (c:y1s) (d:y2s) (i,j) | isGap c || isGap d = diffSSP x1s x2s y1s y2s (i,j)
--same
--diffSSP (x:xs) (y:ys) (i,j) | trace ("Not gap" ++ (show x) ++ " , " ++  (show y) ++ " " ++ (show (i,j))) False = undefined
diffSSP (a:x1s) (b:x2s) (c:y1s) (d:y2s) (i,j) | a==c && b==d = i `seq` j `seq` diffSSP x1s x2s y1s y2s (i+1,j) --same
                                              | a==c = i `seq` j `seq` diffSSP x1s x2s y1s y2s (i+2,j+2) --different 
                                              | a<c = i `seq` j `seq` diffSSP x1s x2s (c:y1s) (d:y2s) (i+1,j+1) --different 
                                              | otherwise  = i `seq` j `seq` diffSSP (a:x1s) (b:x2s) y1s y2s (i+1,j+1) --different  a>c
diffSSP [] [] [] [] (i,j) = (i,j)
diffSSP (a:x1s) (b:x2s) [] [] (i,j) = diffSSP x1s x2s [] [] (i+1,j+1)
diffSSP [] [] (c:y1s) (d:y2s) (i,j) = diffSSP [] [] y1s y2s (i+1,j+1)
