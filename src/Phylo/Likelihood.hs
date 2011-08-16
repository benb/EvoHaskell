module Phylo.Likelihood where
import Phylo.Alignment
import Phylo.Tree
import Phylo.Matrix
import Phylo.Opt
import Numeric.LinearAlgebra ((<>),(<.>),scale,mul,constant,diag)
import Data.Packed.Matrix
import Data.Packed.Vector
import Data.List
import Data.Maybe
import Statistics.Math
import Numeric.GSL.Distribution.Continuous
import Numeric.GSL.Minimization
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
optThmm numCat aln tree pi s = boundedFunction (-1E20) [Nothing,Nothing,Just 0.0] [Nothing,Nothing,Just 1.0] $ partialOptThmm numCat patcounts pi s dataTree where
                                        dataTree = structData numCat AminoAcid pAln transpats tree                                                    
                                        transpats = transpose $ patterns pAln
                                        patcounts = counts pAln
                                        pAln = pAlignment aln

opt2Thmm numCat aln tree pi s = (modelF,dataTree) where
                                dataTree = structData numCat AminoAcid pAln transpats tree                                                                                                                                            
                                transpats = transpose $ patterns pAln  
                                patcounts = counts pAln    
                                pAln = pAlignment aln 
                                modelF = partialOptThmm numCat patcounts pi s
                                

partialOptThmm numCat patcounts pi s dataTree [logalpha,logsigma,priorZero] | trace "OK" True = ans2 where
                                                                                eigenS  = eigQ qMat fullPi                                     
                                                                                ans = logLikelihood patcounts dataTree (flatFullPi numCat pi) eigenS
                                                                                ans2 = ans -- | trace ((show (exp logalpha)) ++ " " ++ (show (exp logsigma)) ++ " " ++ (show priorZero) ++ " " ++ (show ans)) True = ans
                                                                                qMat = thmm numCat pi s priors (exp logalpha) (exp logsigma)
                                                                                remainder=1.0-priorZero
                                                                                priors = priorZero : (replicate numCat (remainder/ (fromIntegral numCat)))
                                                                                fullPi = Data.Packed.Vector.join [mapVector (*priorZero) pi, mapVector (*remainder) $ flatFullPi (numCat-1) pi]


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

getBL node = getBL' node []

getBL' (DLeaf _ bl _ _) bls = bl:bls
getBL' (DINode l r bl) bls = bl:(getBL' l (getBL' r bls))
getBL' (DTree l r) bls = (getBL' l (getBL' r bls))

setBL bls node = fst $ setBL' bls node
setBL' bls (DTree l r) = ((DTree left right), remainder2) where
                        (left,remainder) = setBL' bls l
                        (right,remainder2) = setBL' remainder r

setBL' (bl2:bls) (DINode l r bl) = ((DINode left right bl2), remainder2) where
                             (left,remainder) = setBL' bls l
                             (right,remainder2) = setBL' remainder r

setBL' (bl:bls) (DLeaf a _ b c) = ((DLeaf a bl b c),bls)

optBL:: DNode -> (DNode -> [Double] -> Double) -> [Double] -> [Double] -> Double
optBL dataTree model params bls = model (setBL bls dataTree) params

zeroQMat size = diag $ constant 0.0 size

thmm numCat pi s priors alpha sigma = fixDiag $ fromBlocks subMats where
                                        qMats = (zeroQMat (rows s)) : (map (\(i,mat) -> setRate i mat pi) $ zip (gamma (numCat -1) alpha) $ replicate numCat $ makeQ s pi)
                                        subMats = map (\(i,mat) -> getRow i mat) $ zip [0..] qMats
                                        getRow i mat = map (kk' mat i) [0..(numCat-1)] 
                                        size = cols s
                                        kk' mat i j | i==j = mat
                                                    | otherwise = diagRect 0.0 (mapVector ((priors !! j) * sigma * ) pi) size size


optParamsAndBL (model,dataTree) params lower upper cutoff = optParamsAndBL' model dataTree params lower upper cutoff (model dataTree params) []  where                                                                                                                     

loggedFunc :: ([Double] -> Double) -> ([Double] -> Double)
loggedFunc f = f2 where
               f2 x = ans2 where
                      ans2 | trace ((show x) ++ " -> " ++ (show ans)) True = ans
                      ans = f x
                          
                                                                                                                                                                                                                                              
maximize method pre maxiter params f size = minimize method pre maxiter params f2 size where
                                               f2 x = -(f x)

optParamsAndBL' :: (DNode -> [Double] -> Double) -> DNode -> [Double] -> [Maybe Double] -> [Maybe Double] -> Double -> Double -> [Matrix Double] -> (DNode,[Double],[Matrix Double])
optParamsAndBL' model dataTree params lower upper cutoff lastIter path = optimal where                                                                                                                                                                         
                                                         optBLFunc = loggedFunc $  boundedFunction (-1E20) (repeat $ Just 0.0) (repeat $ Nothing) $ optBL dataTree model optParams                                                                                                                                              
                                                         (optBLs,_) = maximize NMSimplex2 1E-2 1000 (map (\i->0.1) (getBL dataTree)) optBLFunc (getBL dataTree)                                                                                                 
                                                         optDataTree = setBL optBLs dataTree                                                                                                                                               
                                                         model2 = loggedFunc $ boundedFunction (-1E20) lower upper (model dataTree)
                                                         (optParams,_) = maximize NMSimplex2 1E-2 1000 (map (\i->0.1) params) model2 params                                                                        
                                                         thisIter = model optDataTree optParams                                                                                                                                               
                                                         optimal | trace ("Update " ++ (show lastIter) ++ " -> " ++ (show thisIter)) True = if (thisIter-lastIter > cutoff)                                                                                                                                            
                                                                   then optParamsAndBL' model optDataTree optParams lower upper cutoff thisIter path
                                                                   else (optDataTree,optParams,path)   
