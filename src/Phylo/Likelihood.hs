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


-- |Combine two sets of partial likelihoods along two branches
combinePartial :: Matrix Double -> Matrix Double -> Matrix Double
combinePartial = mul 

-- |Combine two lists of partial likelihoods (for mixture models)
allCombine :: [Matrix Double] -> [Matrix Double] -> [Matrix Double]
allCombine (left:lefts) (right:rights) = (combinePartial left right):(allCombine lefts rights)
allCombine [] [] = []

rawPartialLikelihoodCalc :: Matrix Double -> Matrix Double -> Matrix Double
rawPartialLikelihoodCalc pT start = pT <> start

allPLC :: [Matrix Double] -> [Matrix Double] -> [Matrix Double]
allPLC (pT:pTs) (start:starts) = (rawPartialLikelihoodCalc pT start):(allPLC pTs starts)
allPLC [] [] = []


sumLikelihoods :: [Int] -> [Double] -> Double
sumLikelihoods counts likelihoods = foldr foldF 0 $ zip (map fromIntegral counts) likelihoods where
                                        foldF :: (Double,Double) -> Double -> Double 
                                        foldF (i,y) x = x + (i * (log y))
                                        

logLikelihoodModel :: DNode -> Double
logLikelihoodModel (DTree left right pLs pC priors pis) = sumLikelihoods pC likelihoods where
                                                                likelihoodss = map toList $ map (\(pi,pL) -> pi <> pL) $ zip pis pLs
                                                                likelihoods = map summarise (transpose likelihoodss) 
                                                                summarise lkl = sum $ zipWith (*) lkl priors
                                                                
logLikelihood = logLikelihoodModel

data PatternAlignment = PatternAlignment {names :: [String], seqs::[String], columns::[String], patterns::[String], counts::[Int]}

pAlignment (ListAlignment names seqs columns) = PatternAlignment names seqs columns patterns counts where
                                                   (patterns,counts) = unzip $ map (\x -> (head x,length x)) $ group $ sort columns

calcPL :: DNode -> DNode -> [Double] -> [BranchModel] -> [Matrix Double]
calcPL left right dist model = allPLC (map (\x-> x dist) model) $ allCombine (getPL left) (getPL right) 
calcLeafPL :: Matrix Double -> [Double] -> [BranchModel] -> [Matrix Double]
calcLeafPL tips dist model = map (\pT -> rawPartialLikelihoodCalc pT tips) $ (map (\x -> x dist) model)
calcRootPL left right = allCombine (getPL left) (getPL right)


restructData (DLeaf name dist sequence partial _ _) model priors pi  = DLeaf name dist sequence partial model partial' where
                                                                               partial' = calcLeafPL partial [dist] model
restructData (DINode left right dist _ _ ) model priors pi= DINode newleft newright dist model partial where
                                                            partial = calcPL newleft newright [dist] model
                                                            newleft = restructData left model priors pi
                                                            newright = restructData right model priors pi
restructData (DTree left right _ pc _ _ ) model priors pi = DTree newleft newright partial pc priors pi where 
                                            partial = calcRootPL newleft newright
                                            newleft = restructData left model priors pi
                                            newright = restructData right model priors pi

structDataN :: Int -> SeqDataType -> PatternAlignment -> Node -> NNode

structDataN hiddenClasses seqDataType pAln node | trace "OKX" True = structDataN' hiddenClasses seqDataType pAln (transpose $ patterns pAln) node 
structDataN' hiddenClasses seqDataType (PatternAlignment names seqs columns patterns counts) transPat (Leaf name dist) = ans where
                                                                                                                          partial = getPartial hiddenClasses seqDataType sequence 
                                                                                                                          sequence = snd $ fromJust $ find (\(n,_) -> n==name) $ zip names transPat 
                                                                                                                          ans = NLeaf name dist sequence partial 


structDataN' hiddenClasses seqDataType pA transPat (INode c1 c2 dist) = NINode left right dist where
                                                                              left = structDataN' hiddenClasses seqDataType pA transPat c1
                                                                              right = structDataN' hiddenClasses seqDataType pA transPat c2

structDataN' hiddenClasses seqDataType pA transPat (Tree c1 c2 ) = NTree left right $ counts pA where
                                                                              left = structDataN' hiddenClasses seqDataType pA transPat c1
                                                                              right = structDataN' hiddenClasses seqDataType pA transPat c2

structData hiddenClasses seqDataType pAln model transPat node = addModel (structDataN hiddenClasses seqDataType pAln node) model 


addModelNNode :: NNode -> [BranchModel] -> [Double] -> [Vector Double] -> DNode
addModelNNode (NLeaf name dist sequence partial) model _ _  = DLeaf name dist sequence partial model partial' where
                                                                                              partial' = calcLeafPL partial [dist] model

addModelNNode (NINode c1 c2 dist) model priors pi =  DINode left right dist model myPL where
                                                  left = addModelNNode c1 model priors pi
                                                  right = addModelNNode c2 model priors pi
                                                  myPL = calcPL left right [dist] model 
                                                                            
addModelNNode (NTree c1 c2 pat) model priors pi = DTree left right myPL pat priors pi where
                                  left = addModelNNode c1 model priors pi
                                  right= addModelNNode c2 model priors pi
                                  myPL = calcRootPL left right

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

                                 


quickGamma' numCat tree priors pi s aln alpha = logLikelihood dataTree  where
                                                dataTree = structData 1 AminoAcid pAln (map (\r->scaledpT r eigenS) rates) transpats tree priors (repeat pi)
                                                rates = gamma numCat alpha
                                                transpats = transpose $ patterns pAln
                                                eigenS = quickEigen pi s
                                                pAln = pAlignment aln
                                                patcounts = counts pAln
                                                
                                                        
quickEigen pi s = eigQ (normQ (makeQ s pi) pi) pi

standardpT = scaledpT 1.0
scaledpT scale eigenS params = pT eigenS $ (head params) * scale


optGammaF :: Int -> ListAlignment -> Node -> Vector Double -> Matrix Double -> (Double -> Double)
optGammaF numCat aln tree pi s = quickGamma' numCat tree priors pi s aln where
                                                priors = take numCat $ repeat (1.0/(fromIntegral numCat))


--flatFullPi numCat pi = Data.Packed.Vector.mapVector (/(fromIntegral numCat)) $ Data.Packed.Vector.join $replicate numCat pi
flatFullPi numCat pi = makeFullPi (replicate numCat (1.0/(fromIntegral numCat))) pi
makeFullPi priors pi = Data.Packed.Vector.join $ map (\p->Data.Packed.Vector.mapVector (*p) pi) priors
flatPriors numCat = replicate numCat (1.0/ (fromIntegral numCat))

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


type BranchModel = [Double] -> Matrix Double
type OptModel = DNode -> [Double] -> Double
data DNode = DLeaf {dName :: String,dDistance :: Double,sequence::String,tipLikelihoods::Matrix Double,p::([[Double] -> Matrix Double]),partialLikelihoods::[Matrix Double]} |
        DINode DNode DNode Double ([[Double] -> Matrix Double]) [Matrix Double] | 
        DTree DNode DNode [Matrix Double] [Int] [Double] [(Vector Double)]

data NNode = NLeaf String Double String (Matrix Double) | NINode NNode NNode Double | NTree NNode NNode [Int]
data DataModel = DataModel {dTree::DNode, patterncounts::[Int], priors::[Double], pis::[Vector Double]}

class AddTree a where
        addModel :: (AddTree a) => a -> [BranchModel] -> [Double] -> [Vector Double] -> DNode
        addModelFx :: (AddTree a) => a -> ([BranchModel],[Vector Double]) -> [Double] -> DNode
        addModelFx t (bm,pi) prior | trace "-1" True = addModel t bm prior pi
instance AddTree NNode where
        addModel w x y z = ans2 where
                       ans = addModelNNode w x y z
                       ans2 = seq (getPL ans) ans
instance AddTree DNode where
        addModel w x y z = ans2 where 
                           ans = restructData w x y z
                           ans2 = seq (getPL ans) ans

getPL (DLeaf _ _ _ _ _ lkl) = lkl
getPL (DINode _ _ _ _ lkl) = lkl
getPL (DTree _ _ lkl _ _ _) = lkl

getBL node = getBL' node []

getBL' (DLeaf _ bl _ _ _ _) bls = bl:bls
getBL' (DINode l r bl _ _) bls = bl:(getBL' l (getBL' r bls))
getBL' (DTree l r _ _ _ _ ) bls = (getBL' l (getBL' r bls))

setBL bls node = fst $ setBL' bls node
setBL' bls (DTree l r pl pC priors pi) = ((DTree left right partial pC priors pi), remainder2) where
                        (left,remainder) = setBL' bls l
                        (right,remainder2) = setBL' remainder r
                        partial = calcRootPL left right

setBL' (bl2:bls) (DINode l r bl mats pl) = ((DINode left right bl2 mats partial), remainder2) where
                             (left,remainder) = setBL' bls l
                             (right,remainder2) = setBL' remainder r
                             partial = calcPL left right [bl2] mats

setBL' (bl:bls) (DLeaf a _ b tips model _) = ((DLeaf a bl b tips model pl),bls) where
                                      pl = calcLeafPL tips [bl] model 


zeroQMat size = diag $ constant 0.0 size

qdLkl numCat priors dataType aln tree modelF params = logLikelihood $ addModelFx (structDataN numCat dataType (pAlignment aln) tree) (modelF params) priors

type ModelF = [Double] -> ([BranchModel],[Vector Double])
basicModel :: Matrix Double -> Vector Double -> ModelF
basicModel s pi _ = ([model],[pi]) where
                    model = standardpT $ quickEigen pi s

quickLkl aln tree pi s = qdLkl 1 [1.0] AminoAcid aln tree (basicModel s pi) []
quickGamma numCat alpha aln tree pi s = qdLkl 1 (flatPriors numCat) AminoAcid aln tree (gammaModel numCat s pi) [alpha]

gammaModel numCat s pi [alpha] = (models,repeat pi) where
                                models = map (\r -> scaledpT r eigenS) $ gamma numCat alpha
                                eigenS = quickEigen pi s

thmmModel numCat pi s [priorZero,logalpha,logsigma] | traceShow [priorZero,logalpha,logsigma] True = ([thmm numCat pi s priors (exp logalpha) (exp logsigma)],[fullPi]) where
                        priors = priorZero : (replicate (numCat -1) ((1.0-priorZero)/(fromIntegral numCat)) )
                        fullPi = makeFullPi priors pi
 

quickThmm numCat aln tree pi s [priorZero,alpha,sigma] = qdLkl numCat [1.0] AminoAcid aln tree (thmmModel numCat pi s) [priorZero,alpha,sigma] 
                                                        
thmm :: Int -> Vector Double -> Matrix Double -> [Double] -> Double -> Double -> BranchModel
thmm numCat pi s priors alpha sigma = standardpT $ eigQ (fixDiag $ fromBlocks subMats) fullPi where
                                        qMats = (zeroQMat (rows s)) : (map (\(i,mat) -> setRate i mat pi) $ zip (gamma (numCat -1) alpha) $ replicate numCat $ makeQ s pi)
                                        subMats = map (\(i,mat) -> getRow i mat) $ zip [0..] qMats
                                        getRow i mat = map (kk' mat i) [0..(numCat-1)] 
                                        fullPi = makeFullPi priors pi
                                        size = cols s
                                        kk' mat i j | i==j = mat
                                                    | otherwise = diagRect 0.0 (mapVector ((priors !! j) * sigma * ) pi) size size



loggedFunc :: ([Double] -> Double) -> ([Double] -> Double)
loggedFunc f = f2 where
               f2 x = ans3 where
                      ans3 | trace (show x) True = ans2
                      ans2 | trace ((show x) ++ " -> " ++ (show ans)) True = ans
                      ans = f x
                    
maximize method pre maxiter size f params = minimize method pre maxiter size f2 params where
                                               f2 = (*(-1)) . f

optParamsAndBL model tree params priors lower upper cutoff = optParamsAndBL' model tree params priors lower upper cutoff (logLikelihood (addModelFx tree (model params) priors)) []  where                                                                                                                     

optBL :: (AddTree t) => ModelF -> t -> [Double] -> [Double] -> Double -> DNode
optBL model t' params priors cutoff = optTree where
        tree = addModelFx t' (model params) priors
        optBLFunc = loggedFunc . boundedFunction (-1E20) (repeat $ Just 0.0) (repeat $ Nothing) . rawOptBLFunc   
        rawOptBLFunc t params = logLikelihood $ setBL params t
        (optBLs,_) = maximize NMSimplex2 1E-3 10000 (map (\i->0.1) (getBL tree)) (optBLFunc tree) (getBL tree)
        optTree = setBL optBLs tree
 
optParams :: (AddTree t) => ModelF -> t -> [Double] -> [Double] -> [Maybe Double] -> [Maybe Double] -> Double -> [Double]
optParams model t' params priors lower upper cutoff = bestParams where
        tree = addModelFx t' (model params) priors
        optParamsFunc = loggedFunc . boundedFunction (-1E20) lower upper . rawOptParamsFunc       
        rawOptParamsFunc t params = logLikelihood $ addModelFx t (model params) priors                                                                                    
        (bestParams,_) = maximize NMSimplex2 1E-4 100000 (map (\i->0.1) params) (optParamsFunc tree) params    
        


optParamsAndBL' :: (AddTree t) => ModelF -> t -> [Double] -> [Double] -> [Maybe Double] -> [Maybe Double] -> Double -> Double -> [Matrix Double] -> (DNode,[Double],[Matrix Double])  
optParamsAndBL' model tree params priors lower upper cutoff lastIter path | trace ("OK2" ++ (show lastIter)) $ True = optimal where                                                                                                                                                        
                                              bestParams = optParams model tree params priors lower upper cutoff
                                              bestTree = optBL model tree bestParams priors cutoff
                                              optimal = (addModelFx bestTree (model bestParams) priors,bestParams,path)


                                               {-numBL = length $ getBL dataTree                                                                                                                                                      -}
                                               {-optParams2 = optBLs ++ optParams                                                                                                                                                     -}
                                               {-model3 = loggedFunc $ boundedFunction (-1E20) ((replicate numBL $ Just 0.0) ++ lower) ((replicate numBL $ Nothing) ++ upper) $ optBLBoth dataTree model                              -}
                                               {-(optParams3,_) = maximize NMSimplex2 1E-4 100000 (map (\i->0.1) optParams2) model3 optParams2                                                                                        -}
                                               {-bestBL = optParams3                                                                                                                                                                  -}
                                               {-(bestTree,bestParams)=setBL' bestBL dataTree                                                                                                                                         -}
                                               {-thisIter = model bestTree bestParams                                                                                                                                                 -}
                                               {-optimal | trace ("Update " ++ (show lastIter) ++ " -> " ++ (show thisIter)) True = if (thisIter-lastIter > cutoff)                                                                   -}
                                                         {-then optParamsAndBL' model bestTree bestParams lower upper cutoff thisIter path                                                                                            -}
                                                         {-else (bestTree,bestParams,path)      -}                 
