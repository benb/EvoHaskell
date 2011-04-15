module Phylo.Likelihood where
import Phylo.Alignment
import Phylo.Tree
import Phylo.Matrix
import Numeric.LinearAlgebra ((<>),(<.>))
import Data.Packed.Matrix
import Data.Packed.Vector
import Data.List
import Data.Maybe

partialLikelihood' ::  (Double -> [Vector Double] -> [Vector Double]) -> DNode -> [Vector Double]
partialLikelihood' f (DLeaf name dist sequence partial) = f dist partial
partialLikelihood' f (DINode c1 c2 dist) = f dist $ combinePartial (partialLikelihood' f c1) (partialLikelihood' f c2)
partialLikelihood' f (DTree c1 c2) = combinePartial (partialLikelihood' f c1) (partialLikelihood' f c2)

combinePartial :: [Vector Double] -> [Vector Double] -> [Vector Double]
combinePartial a b = zipWith (zipVectorWith (*)) a b

partialLikelihoodCalc :: EigenS -> Double -> [Vector Double] -> [Vector Double]
partialLikelihoodCalc eig t pL = map (myPt <>) pL where
                                        myPt = pT eig t

partialLikelihood eigenS = partialLikelihood' $ partialLikelihoodCalc eigenS
                                        
logLikelihood :: [Int] -> DNode -> Vector Double -> EigenS -> Double
logLikelihood counts dataTree pi eigenS = foldr foldF 0 $ zip (map fromIntegral counts) likelihoods where
                                             likelihoods = map (pi <.>) pL
                                             pL = partialLikelihood eigenS dataTree
                                             foldF :: (Double,Double) -> Double -> Double 
                                             foldF (i,y) x = x + (i * (log y))

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
getPartial :: Int -> SeqDataType -> String -> [Vector Double]
getPartial _ _ [] = []
getPartial classes AminoAcid (x:xs) = (aaPartial classes x):(getPartial classes AminoAcid xs)
getPartial classes Nucleotide (x:xs) = error "Nucleotides unimplemented"

aaOrder = ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
aaPartial classes x | isGapChar x = buildVector (20*classes) (\i->1.0) 
                    | otherwise = case findIndex (==x) aaOrder of 
                                       Just index -> buildVector (20*classes) (\i -> if i`mod`20==index then 1.0 else 0.0) where
                                       Nothing -> error $ "character " ++ [x] ++ "is not an amino acid or gap"

                                 
quickLkl aln tree pi s = logLikelihood patcounts dTree pi eigenS where
                                dTree = structData 1 AminoAcid pAln transpats tree
                                pAln = pAlignment aln
                                transpats = transpose $ patterns pAln
                                qMat = normQ (makeQ s pi) pi 
                                eigenS = eigQ qMat pi
                                patcounts = counts pAln

