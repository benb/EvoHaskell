module Alignment.Test where
import Test.QuickCheck
import Data.List
import Alignment
import qualified Alignment.Dist as DF
import qualified Alignment.DistSet as DS
import Control.Monad
import Debug.Trace
import Text.Printf

main  = mapM_ (\(s,a) -> printf "%-25s: " s >> a) alltests 

deepCheck = quickCheckWith (stdArgs { maxSuccess = 1000})


alltests = [("sort",deepCheck prop_sortaln),
            ("nonneg",quickCheck prop_nonneg),
            ("discriminant",quickCheck prop_disc),
            ("symmetric",quickCheck prop_sym),
            ("triangle ineq",quickCheck prop_trig),
            ("same as sets",quickCheck prop_sets)]


toDist :: (Int,Int) -> Double
toDist (i,j) = (fromIntegral j)/(fromIntegral i)

allDists = [DF.homDist,DF.hom0Dist,DF.homGapDist]
allSetDists = [DS.homDist, DS.hom0Dist, DS.homGapDist]

test :: [ListAlignment -> ListAlignment -> (Int,Int)] -> ListAlignment -> ListAlignment -> ((Int,Int) -> Bool) -> Bool
test testSet x y criteria = or $ map (\f-> criteria(f x y)) testSet


test2 :: [ListAlignment -> ListAlignment -> (Int,Int)] -> ListAlignment -> ListAlignment -> ListAlignment -> ListAlignment -> ((Int,Int) -> (Int,Int) -> Bool) -> Bool
test2 testSet x y x' y' criteria = or $ map (\f-> criteria (f x y) (f x' y')) testSet

testTriple :: [ListAlignment -> ListAlignment -> (Int,Int)] -> ListAlignment -> ListAlignment -> ListAlignment -> ((Int,Int) -> (Int,Int) -> (Int,Int) -> Bool) -> Bool
testTriple testSet x y z criteria = or $ map (\f-> criteria (f x y) (f x z) (f y z)) testSet

testZip :: [ListAlignment -> ListAlignment -> (Int,Int)] -> [ListAlignment -> ListAlignment -> (Int,Int)] -> ListAlignment -> ListAlignment -> ((Int,Int) -> (Int,Int) -> Bool) -> Bool
testZip first second x y criteria = or $ map (\(i,j)-> criteria i j) (firstAns `zip` secondAns) where
                                         firstAns = map (\f-> f x y) first
                                         secondAns = map (\f-> f x y) second

prop_nonneg (AlignmentPair x y) = test allDists x y (\x-> toDist x >=0) 
prop_disc (AlignmentPair x y) = test allDists x y (\z-> (x==y) || (toDist z)>0) 
prop_sym (AlignmentPair x y) = test2 allDists x y y x (\a b-> a==b)
prop_trig (AlignmentTriple x y z) = testTriple allDists x y z (\a b c -> ((toDist a) + (toDist b) >= (toDist c))) -- (toDist $ DF.homDist x y) + (toDist $ DF.homDist y z) >= (toDist $ DF.homDist x z)

prop_sets (AlignmentPair x y) = testZip allDists allSetDists x y (\x y-> x==y)

prop_sortaln (AlignmentPair x y) = (map (filter (not.isGapChar)) (sequences (sortAlignment x))) == (map (filter (not.isGapChar)) (sequences x))


data UngappedSequence = UngappedSequence [SeqChar] deriving Show


getSeq (UngappedSequence x) = map getSeqChar x

data SeqChar = SeqChar Char deriving Show
getSeqChar (SeqChar x) = x
instance Arbitrary SeqChar where
        arbitrary = liftM SeqChar $ elements ['A'..'Z']

instance Arbitrary UngappedSequence where
        arbitrary = do x <- listOf1 (arbitrary :: Gen SeqChar)
                       return $ UngappedSequence x

data AlignmentPair = AlignmentPair ListAlignment ListAlignment deriving Show
data AlignmentTriple = AlignmentTriple ListAlignment ListAlignment ListAlignment deriving Show
       
instance Arbitrary AlignmentPair where
        arbitrary = do unalignedSeqs <- liftM (map getSeq) $ (liftM2 (++)) (listOf1 (arbitrary :: Gen UngappedSequence)) (listOf1 (arbitrary :: Gen UngappedSequence))
                       finalLen <- do x<-(arbitrary :: Gen Int)
                                      let ungaplen = (maximum $ map length unalignedSeqs)
                                      return $ ((abs x) `mod` ungaplen) + ungaplen 
                       shuffled1 <- liftM (reapplySeq unalignedSeqs) $ mapM shuffle $ pad unalignedSeqs finalLen
                       shuffled2 <- liftM (reapplySeq unalignedSeqs) $ mapM shuffle $ pad unalignedSeqs finalLen
                       return $ AlignmentPair (safeListAlignment (names shuffled1) shuffled1) (safeListAlignment (names shuffled2) shuffled2) where
                       pad [] i = []
                       pad (x:xs) i = (x ++ (take (i-(length x)) $ repeat '-')):(pad xs i)
                       names = map (\x -> [fst x]) . zip ['A'..] 

reapplySeq (x:xs) (y:ys) = (reapplySeq' x y) : reapplySeq xs ys
reapplySeq [] [] = []
reapplySeq' (x:xs) (y:ys) 
           | isGapChar y = y : reapplySeq' (x:xs) ys
           | otherwise = x : reapplySeq' xs ys 
reapplySeq' [] (y:ys) = y:ys
reapplySeq' [] [] = [] 

shuffle :: (Eq a) => [a] -> Gen [a]
shuffle [] = return []
shuffle xs = do x  <- oneof $ map return xs
                ys <- shuffle $ delete x xs
                return (x:ys)

instance Arbitrary AlignmentTriple where
        arbitrary = do unalignedSeqs <- liftM (map getSeq) $ (liftM2 (++)) (listOf1 (arbitrary :: Gen UngappedSequence)) (listOf1 (arbitrary :: Gen UngappedSequence))
                       finalLen <- do x<-(arbitrary :: Gen Int)
                                      return $ ((abs x) `mod` 10) + (maximum $ map length unalignedSeqs)
                       shuffled1 <- liftM (reapplySeq unalignedSeqs) $ mapM shuffle $ pad unalignedSeqs finalLen
                       shuffled2 <- liftM (reapplySeq unalignedSeqs) $ mapM shuffle $ pad unalignedSeqs finalLen
                       shuffled3 <- liftM (reapplySeq unalignedSeqs) $ mapM shuffle $ pad unalignedSeqs finalLen
                       return $ AlignmentTriple (safeListAlignment (names shuffled1) shuffled1) (safeListAlignment (names shuffled2) shuffled2) (safeListAlignment (names shuffled3) shuffled3) where
                       pad [] i = []
                       pad (x:xs) i = (x ++ (take (i-(length x)) $ repeat '-')):(pad xs i)
                       names = map (\x -> [fst x]) . zip ['A'..] 
