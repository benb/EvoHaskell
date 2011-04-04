module Alignment.DistSet where
import Alignment
import Debug.Trace
import qualified Data.Set as Set


hom0Dist = zeroDist numberifyBasic
homDist = labDist (setPair []) setDist numberifyBasic
homGapDist = labDist (setPair []) setDist numberifyGap
homTreeDist t = labDist (setPair []) setDist (numberifyGapTree t)

isPermutation :: ListAlignment -> ListAlignment -> Bool
isPermutation a b | (names a) /= (names b) = False
isPermutation a b = (dropGaps a) == (dropGaps b)


genDist :: (ListAlignment -> [[(Int)]]) -> ListAlignment -> ListAlignment -> (Int,Int)
genDist f = labDist (setPair []) setDist f

zeroDist f = labDist (setPairTrig) setDistZero f
                        
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



-- |'labDist' computes the distance between two alignments after labelling
-- it takes a labelling function and two alignments
-- and returns a a tuple of (denonimator,numerator), i.e. distance is
-- snd/fst
labDist :: (Ord a,SiteLabel a) => ([[(a,Int)]] -> [[[((a,Int),(a,Int))]]]) -> (Set.Set ((a,Int),(a,Int))-> Set.Set ((a,Int),(a,Int)) -> (Int,Int)) -> (ListAlignment -> [[(a)]]) -> ListAlignment -> ListAlignment -> (Int,Int)
labDist pairFunc distFunc numF aln1 aln2 | trace "Set Distances " False = undefined
labDist pairFunc distFunc numF aln1 aln2 = distFunc set1 set2 where
                           setLabel :: [[a]] -> Int -> [[(a,Int)]]
                           setLabel [] i = []
                           setLabel (x:xs) i = (map (\f -> (f,i)) x) :  (setLabel xs (i+1))
                           set1 = (toSet (pairFunc (setLabel (numF aln1) 0)))
                           set2 = (toSet (pairFunc (setLabel (numF aln2) 0)))
                           toSet (x:xs) = (toSet2 x) `Set.union` (toSet xs)
                           toSet [] = Set.fromList []
                           toSet1 (x:xs) = (toSet2 x) `Set.union` (toSet1 xs)
                           toSet1 [] = Set.fromList []
                           toSet2 (x:xs) = (Set.fromList x) `Set.union` (toSet2 xs)
                           toSet2 [] = Set.fromList []

setPairTrig [] = []
setPairTrig (x:xs) = (map (\l -> filter neitherGap $ l `zip` x) xs ) : setPairTrig xs

--need to filter out (gap,x) comparisons

isNotGap :: SiteLabel a => ((a,b),(a,b)) -> Bool
isNotGap ((i,j),x) = not $ isGap i

neitherGap :: SiteLabel a => ((a,b),(a,b)) -> Bool
neitherGap ((i,j),(k,l)) = (not ((isGap i) || (isGap k)))

setPair xs [] =  []
setPair ys (x:xs) = (map (\l -> filter isNotGap $ l `zip` x) (xs ++ ys)) : setPair (x:ys) xs


setDist :: (SiteLabel a, Ord a) => Set.Set ((a,Int),(a,Int))-> Set.Set ((a,Int),(a,Int)) -> (Int,Int)
--setDist a b | trace (show a) False = undefined
--setDist a b | trace (show b) False = undefined
setDist a b = ((Set.size a)+(Set.size b),Set.size ((a Set.\\ b) `Set.union` (b Set.\\ a)))

setDistZero a b = ((Set.size $ a `Set.union` b), Set.size ((a Set.\\ b) `Set.union` (b Set.\\ a)))
