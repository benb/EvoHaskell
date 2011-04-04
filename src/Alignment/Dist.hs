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
labDistSeq :: ([(a,a)]->[(a,a)]->(Int,Int)->(Int,Int))-> [a] -> [a] -> [[a]] ->  [[a]] -> (Int,Int) -> (Int,Int)
labDistSeq f seqA seqB (seqA2:xs) (seqB2:ys) ans = labDistSeq f seqA seqB xs ys $! (f (zip seqA seqA2) (zip seqB seqB2) ans)
labDistSeq f seqA seqB [] [] ans = ans
       
-- | compute distance for pairs of labels 
diffIn :: (SiteLabel a)=> [(a,a)] -> [(a,a)] -> (Int,Int) -> (Int,Int)
--diffIn (x:xs) (y:ys) ans | trace ((show x) ++ (show y) ++ (show ans)) False = undefined
--First, skip gaps on left side of xs and ys
diffIn ((x1,x2):xs) y (i,j) | (isGap x1) = i `seq` j `seq` diffIn xs y (i,j)
diffIn x ((y1,y2):ys) (i,j) | (isGap y1) = i `seq` j `seq` diffIn x ys (i,j)

--Metric 0, ignore gaps
--diffIn ((x1,(x2i,x2g)):xs) ((y1,(y2i,y2g)):ys) (i,j) | (x2i<0 && y2i>0 && x2g==Nothing) || (x2i>0 && y2i<0 && y2g==Nothing)  = i `seq` j `seq` diffIn xs ys  (i+1,j+1)
--Same
diffIn ((x1,x2):xs) ((y1,y2):ys) (i,j) | x2==y2  = i `seq` j `seq` diffIn xs ys  (i+2,j)
--Different
diffIn ((x1,x2):xs) ((y1,y2):ys) (i,j)  = i `seq` j `seq` diffIn xs ys (i+2,j+2) 
--diffIn xs ys ans | trace ((show xs) ++ (show ys)) False = undefined
diffIn [] [] t = t



-- | compute distance for pairs of labels for metric 0 (SSP)
diffSSP :: (SiteLabel a, Ord a) => [(a,a)] -> [(a,a)] -> (Int,Int) -> (Int,Int)
--diffSSP (x:xs) (y:ys) (i,j) | trace ("diffSSP" ++ (show x) ++ " , " ++  (show y) ++ " " ++ (show (i,j))) False = undefined
diffSSP ((a,b):xs) (ys) (i,j) | isGap a || isGap b = diffSSP xs ys (i,j)
diffSSP (xs) ((c,d):ys) (i,j) | isGap c || isGap d = diffSSP xs ys (i,j)
--same
--diffSSP (x:xs) (y:ys) (i,j) | trace ("Not gap" ++ (show x) ++ " , " ++  (show y) ++ " " ++ (show (i,j))) False = undefined
diffSSP ((a,b):xs) ((c,d):ys) (i,j) | a==c && b==d = i `seq` j `seq` diffSSP xs ys (i+1,j) --same
--diffSSP (x:xs) (y:ys) (i,j) | trace ("Not same" ++ (show x) ++ " , " ++  (show y)) False = undefined
diffSSP ((a,b):xs) ((c,d):ys) (i,j) | a==c = i `seq` j `seq` diffSSP xs ys (i+2,j+2) --different 
diffSSP ((a,b):xs) ((c,d):ys) (i,j) | a<c = i `seq` j `seq` diffSSP xs ((c,d):ys) (i+1,j+1) --different 
diffSSP ((a,b):xs) ((c,d):ys) (i,j)  = i `seq` j `seq` diffSSP ((a,b):xs) ys (i+1,j+1) --different  a>c
--diffSSP [] [] (i,j) | trace "Ping" False  = undefined
diffSSP [] [] (i,j) = (i,j)
diffSSP (x:xs) [] (i,j) = diffSSP xs [] (i+1,j+1)
diffSSP [] (y:ys) (i,j) = diffSSP [] ys (i+1,j+1)

