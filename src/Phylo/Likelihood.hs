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
import Debug.Trace


partialLikelihood' ::  (DNode -> Matrix Double -> Matrix Double) -> DNode -> Matrix Double
partialLikelihood' f (DLeaf name dist sequence partial) = f (DLeaf name dist sequence partial) partial
partialLikelihood' f (DINode c1 c2 dist) = par p1 $ f (DINode c1 c2 dist) $ combinePartial p1 p2 where
                                                p1 = partialLikelihood' f c1
                                                p2 = partialLikelihood' f c2
partialLikelihood' f (DTree c1 c2) = par p1 $ combinePartial p1 p2 where
                                                p1 = partialLikelihood' f c1
                                                p2 = partialLikelihood' f c2

combinePartial :: Matrix Double -> Matrix Double -> Matrix Double
combinePartial = mul 

--this is for time-homogenous models
partialLikelihoodCalc :: EigenS -> DNode -> Matrix Double -> Matrix Double
partialLikelihoodCalc eig (DLeaf _ t _ _) = partialLikelihoodCalc' eig t
partialLikelihoodCalc eig (DINode _ _ t) = partialLikelihoodCalc' eig t
partialLikelihoodCalc eig (DTree _ _) = partialLikelihoodCalc' eig 0.0
partialLikelihoodCalc' eig t mypL = myPt <> mypL where
                                        myPt = pT eig t

partialLikelihoodCalcHet :: (DNode -> EigenS) -> DNode -> Matrix Double -> Matrix Double
partialLikelihoodCalcHet map node = case node of 
                                        (DLeaf _ t _ _) -> partialLikelihoodCalc' (map node) t
                                        (DINode _ _ t)  -> partialLikelihoodCalc' (map node) t
                                        (DTree _ _)     -> partialLikelihoodCalc' (map node) 0.0
                                        
 
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


flatFullPi numCat pi = Data.Packed.Vector.mapVector (/(fromIntegral numCat)) $ Data.Packed.Vector.join $replicate numCat pi
optGammaF2 :: Int -> ListAlignment -> Node -> Vector Double -> Matrix Double -> Double -> Double
optGammaF2 numCat aln tree pi s alpha = logLikelihood patcounts dataTree fullPi eigenS where
                                                dataTree = structData numCat AminoAcid pAln transpats tree
                                                fullPi = flatFullPi numCat pi
                                                pAln = pAlignment aln
                                                transpats = transpose $ patterns pAln
                                                patcounts = counts pAln
                                                qMat = normQ (makeQ s pi) pi
                                                rateMats = map (\(i,mat) -> setRate i mat pi ) $ zip rates $ replicate numCat qMat
                                                rates = gamma numCat alpha
                                                fullMat = combineQ rateMats
                                                eigenS = eigQ fullMat fullPi

optThmm :: Int -> ListAlignment -> Node -> Vector Double -> Matrix Double -> ([Double] -> Double)
optThmm numCat aln tree pi s = partialOptThmm numCat patcounts dataTree pi s where
                                        dataTree = structData numCat AminoAcid pAln transpats tree                                                    
                                        transpats = transpose $ patterns pAln
                                        patcounts = counts pAln
                                        pAln = pAlignment aln


partialOptThmm numCat patcounts dataTree pi s [logalpha,logsigma] = ans2 where
                                                                                eigenS  = eigQ qMat fullPi                                     
                                                                                ans = logLikelihood patcounts dataTree (flatFullPi numCat pi) eigenS
                                                                                ans2 | trace ((show (exp logalpha)) ++ " " ++ (show (exp logsigma)) ++ " " ++ (show ans)) True = ans
                                                                                qMat = thmm numCat pi s priors (exp logalpha) (exp logsigma)
                                                                                priors = replicate numCat (1.0 / (fromIntegral numCat))
                                                                                fullPi = flatFullPi numCat pi

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


data DNode = DLeaf {dName :: String,dDistance :: Double,sequence::String,likelihoods::Matrix Double} | DINode DNode DNode Double | DTree DNode DNode 

thmm numCat pi s priors alpha sigma = fixDiag $ fromBlocks subMats where
                                        qMats = map (\(i,mat) -> setRate i mat pi) $ zip (gamma numCat alpha) $ replicate numCat $ makeQ s pi
                                        subMats = map (\(i,mat) -> getRow i mat) $ zip [0..] qMats
                                        getRow i mat = map (kk' mat i) [0..(numCat-1)] 
                                        size = cols s
                                        kk' mat i j | i==j = mat
                                                    | otherwise = diagRect 0.0 (mapVector ((priors !! j) * sigma * ) pi) size size
