module Phylo.Likelihood where
import Phylo.Alignment
import Phylo.Tree
import Phylo.Matrix
import Numeric.LinearAlgebra ((<>),(<.>),scale,mul)
import Data.Packed.Matrix
import Data.Packed.Vector
import Data.List
import Data.Maybe
import Statistics.Math
import Numeric.GSL.Distribution.Continuous
import Control.Exception as E
import Control.Parallel
import Control.Parallel.Strategies


partialLikelihood' ::  (Double -> Matrix Double -> Matrix Double) -> DNode -> Matrix Double
partialLikelihood' f (DLeaf name dist sequence partial) = f dist partial
partialLikelihood' f (DINode c1 c2 dist) = par p1 $ f dist $ combinePartial p1 p2 where
                                                p1 = partialLikelihood' f c1
                                                p2 = partialLikelihood' f c2
partialLikelihood' f (DTree c1 c2) = par p1 $ combinePartial p1 p2 where
                                                p1 = partialLikelihood' f c1
                                                p2 = partialLikelihood' f c2

combinePartial :: Matrix Double -> Matrix Double -> Matrix Double
combinePartial = mul 

partialLikelihoodCalc :: EigenS -> Double -> Matrix Double -> Matrix Double
partialLikelihoodCalc eig t mypL = myPt <> mypL where
                                        myPt = pT eig t
 
partialLikelihood :: EigenS ->  DNode -> Matrix Double
partialLikelihood eigenS = partialLikelihood' $ partialLikelihoodCalc eigenS
                                        
logLikelihood :: [Int] -> DNode -> Vector Double -> EigenS -> Double
logLikelihood counts dataTree pi eigenS = sumLikelihoods counts likelihoods where
                                               likelihoods = toList $ pi <> pL
                                               pL = partialLikelihood eigenS dataTree

sumLikelihoods :: [Int] -> [Double] -> Double
sumLikelihoods counts likelihoods = foldr foldF 0 $ zip (map fromIntegral counts) likelihoods where
                                        foldF :: (Double,Double) -> Double -> Double 
                                        foldF (i,y) x = x + (i * (log y))
                                        

logLikelihoodMixture counts dataTree priors pis eigenS = sumLikelihoods counts likelihoods where
                                                                pLs = (parMap rseq) (\x-> partialLikelihood x dataTree) eigenS
                                                                likelihoodss = (parMap rseq) (\(pi,pL)-> toList $ pi <> pL) $ zip pis pLs
                                                                likelihoods = map summarise (transpose likelihoodss) 
                                                                summarise lkl = sum $ zipWith (*) lkl priors

gammaMix numCat alpha (u,lambda,u') = map (\s -> (u,scale s lambda,u')) scales where
                           scales = gamma numCat alpha

gamma :: Int -> Double -> [Double]
gamma numCat shape | shape > 50000.0 = gamma numCat 50000.0 --work around gsl convergence errors (who wants >100 in the gamma dist anyway?)
                   | otherwise       = map rK' [0..(numCat-1)] where
                        alpha = shape
                        beta = shape
                        factor = alpha/beta*(fromIntegral numCat)

                        freqK = map freqKf [0..(numCat-1)]
                        freqKf i = incompleteGamma (alpha + 1.0) (beta * (gammaInvCDF (((fromIntegral i)+1.0)/(fromIntegral numCat))))

                        rK' 0 = (head freqK) * factor
                        rK' n | n < numCat = ((freqK!!n) - (freqK!!(n-1))) * factor
                              | otherwise = factor * (1.0-(freqK!!(numCat-2)))
                        gammaInvCDF p = (density_1p ChiSq UppInv (2.0*alpha) (1.0-p)) / (2.0*beta) -- (UppInv mysteriously has more stable convergence)


data PatternAlignment = PatternAlignment {names :: [String], seqs::[String], columns::[String], patterns::[String], counts::[Int]}

pAlignment (ListAlignment names seqs columns) = PatternAlignment names seqs columns patterns counts where
                                                   (patterns,counts) = unzip $ map (\x -> (head x,length x)) $ group $ sort columns

structData :: Int -> SeqDataType -> PatternAlignment -> [String] -> Node -> DNode
structData hiddenClasses seqDataType (PatternAlignment names seqs columns patterns counts) transPat (Leaf name dist) = DLeaf name dist sequence partial where
                                                                                                                    sequence = snd $ fromJust $ find (\(n,_) -> n==name) $ zip names transPat 
                                                                                                                    partial = getPartial hiddenClasses seqDataType sequence

structData hiddenClasses seqDataType pA transPat (INode c1 c2 dist) = DINode (structData hiddenClasses seqDataType pA transPat c1) (structData hiddenClasses seqDataType pA transPat c2) dist 
structData hiddenClasses seqDataType (PatternAlignment names seqs columns patterns counts) transPat (Tree c1 c2) = DTree (structData hiddenClasses seqDataType pA transPat c1) (structData hiddenClasses seqDataType pA transPat c2) where
                                                                                                                        pA = PatternAlignment names seqs columns patterns counts

data SeqDataType = AminoAcid | Nucleotide
getPartial :: Int -> SeqDataType -> String -> Matrix Double
getPartial a b c = trans $ fromRows $ getPartial' a b c
getPartial' :: Int -> SeqDataType -> String -> [Vector Double]
getPartial' _ _ [] = []
getPartial' classes AminoAcid (x:xs) = (aaPartial classes x):(getPartial' classes AminoAcid xs)
getPartial' classes Nucleotide (x:xs) = error "Nucleotides unimplemented"

aaOrder = ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
aaPartial classes x | isGapChar x = buildVector (20*classes) (\i->1.0) 
                    | otherwise = case findIndex (==x) aaOrder of 
                                       Just index -> buildVector (20*classes) (\i -> if i`mod`20==index then 1.0 else 0.0) where
                                       Nothing -> error $ "character " ++ [x] ++ "is not an amino acid or gap"

                                 
quickLkl aln tree pi s = logLikelihood patcounts dTree pi eigenS where
                                dTree = structData 1 AminoAcid pAln transpats tree
                                pAln = pAlignment aln
                                transpats = transpose $ patterns pAln
                                eigenS = quickEigen pi s
                                patcounts = counts pAln

quickGamma' numCat patcounts dataTree priors pi eigenS alpha = logLikelihoodMixture patcounts dataTree priors (repeat pi) eigenSs where
                                                        eigenSs = gammaMix numCat alpha eigenS

quickGamma numCat alpha aln tree pi s  = optGammaF numCat aln tree pi s alpha

                                                
quickEigen pi s = eigQ (normQ (makeQ s pi) pi) pi


optGammaF :: Int -> ListAlignment -> Node -> Vector Double -> Matrix Double -> (Double -> Double)
optGammaF numCat aln tree pi s = quickGamma' numCat patcounts dataTree priors pi eigenS where
                                                dataTree = structData 1 AminoAcid pAln transpats tree
                                                pAln = pAlignment aln
                                                transpats = transpose $ patterns pAln
                                                priors = take numCat $ repeat (1.0/(fromIntegral numCat))
                                                eigenS = quickEigen pi s
                                                patcounts = counts pAln

data DNode = DLeaf {dName :: String,dDistance :: Double,sequence::String,likelihoods::Matrix Double} | DINode DNode DNode Double | DTree DNode DNode 
