{-# Language ScopedTypeVariables,CPP #-}
module Phylo.Likelihood (optBSParamsBL,pAlignment,logLikelihood,
       DNode(DLeaf,DINode,DTree),PatternAlignment(PatternAlignment),addModelFx,structDataN
       ,mapBack,toPBELengths,getQMat,rawlikelihoods,getPi,getPriors,toPBEList,getPartialBranchEnds
       ,gammaModel,annotateTreeWith,flatPriors,thmmModel,optBLDFull0,SeqDataType(AminoAcid,Nucleotide),gammaModelQ,thmmModelQ
       ,getLeftSplit,leftSplit,cachedBranchModelTree,makeSimulatedAlignment,getAllF,makeMapping,makeSimulatedAlignmentWithGaps,Phylo.Likelihood.columns,genList
       ,addModelNNode,removeModel,quickThmm,quickGamma,thmmPerBranchModel,annotateTreeWithNumberSwitches,annotateTreeWithNumberSwitchesSigma
       ,jcF,gtrF,gtrS,wagF,jttF,hkyF,hkyS,dataSize,customF,getSensibleParams,SPiFunctionTuple,zeroParam,piByLog) where
import Phylo.Alignment
import Phylo.Tree
import Phylo.Matrix
import Phylo.Opt
import Phylo.NLOpt
import Numeric.LinearAlgebra ((<>),(<.>),scale,mul,constant,diag,kronecker,add,ident)
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
import Phylo.Data
import System.Random
import Numeric.GSL.Differentiation
import Numeric.GSL.OneDimensionalMinimization
import qualified Data.Vector as DVec
import Data.Vector ((!))
import Numeric.LinearAlgebra.LAPACK
import Control.Parallel.Strategies
import Control.DeepSeq
import Data.Ord
import qualified Data.HashMap as HM
import Control.Monad.ST
import Control.Monad
import Data.Packed.ST

#ifdef DEBUG
mytrace = trace
#else
mytrace x y = y
#endif

mytraceShow = mytrace . show

dataSize AminoAcid = 20
dataSize Nucleotide = 4


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
logLikelihoodModel root@(DTree _ _ _ _ pC _ _) = sumLikelihoods pC $ likelihoods root where

concatNewLine :: [[String]] -> String
concatNewLine = intercalate "\n" . map (intercalate "\n")


posteriorTipsCSV :: Int -> PatternAlignment -> DNode -> String
posteriorTipsCSV i aln j = intercalate "\n" $ (posteriorTipsCSVX aln) $ posteriorTips i aln j

posteriorTipsCSVX aln = map (posteriorTipsCSVX' aln)

cachedBranchModel initDist bm = let cached = (bm initDist) in
                                \x -> case x of 
                                           x | x == initDist -> cached
                                             | otherwise -> bm x

cachedBranchModelTree (DTree l m r pLs pC priors pis) = DTree (cachedBranchModelTree l) (cachedBranchModelTree m) (cachedBranchModelTree r) pLs pC priors pis
cachedBranchModelTree (DINode l r bl model pl) = DINode (cachedBranchModelTree l) (cachedBranchModelTree r) bl (map (cachedBranchModel bl) model) pl
cachedBranchModelTree (DLeaf name dist seq partial model pl) = DLeaf name dist seq partial (map (cachedBranchModel dist) model) pl

posteriorTipsCSVX' :: PatternAlignment ->  [(Int,String,DNode,[Double])] -> String
posteriorTipsCSVX' (PatternAlignment names seqs _ _ _)  list = intercalate "\n" $ map (\(x,x',y,z) -> name ++ "," ++  (show x) ++ "," ++ (show x') ++"," ++ [z]  ++"," ++ (intercalate "," (map show y))) $ zip4 [0..] realseqpos (transpose post) alignedseq where
                         post = map (\(x,y,z,lkl)->lkl) list
                         name = head $ map (\(x,y,z,lkl) -> y) list
                         id = fromJust $ findIndex (==name) names
                         alignedseq = seqs !! id
                         realseqpos = reverse $ foldl incPos [] alignedseq
                         incPos [] '-' = [-1]
                         incPos [] base = [0]
                         incPos (x:xs) '-' = (x:x:xs) 
                         incPos (x:xs) base = (x+1:x:xs)


posteriorTipsCSV' :: (Int,String,DNode,[Double]) -> String
posteriorTipsCSV' (i,j,k,lkl) = j ++ "," ++ (show i) ++ "," ++ (intercalate "," $ map show lkl)

posteriorTips :: Int -> PatternAlignment -> DNode -> [[(Int,String,DNode,[Double])]]
posteriorTips numModels aln root@(DTree l m r pLs pC priors pis) = map normaliseL2 $ posteriorTipsLkl numModels aln root where
                                                                --   getLkl (i,j,k,lkl) = lkl
                                                                --   getLklList = (map . map) getLkl
normaliseL2 :: [(Int,String,DNode,[Double])] -> [(Int,String,DNode,[Double])]
normaliseL2 list = replaced where
                    unzipped = unzip4 list
                    (_,_,_,lkls) = unzipped
                    probs = normaliseL2' lkls
                    replaced = map (\((i,j,k,lkl),probs) -> (i,j,k,probs)) $ zip list probs

normaliseL2' :: [[Double]] -> [[Double]]
normaliseL2' list = transpose $ normaliseL2'' $ transpose list
normaliseL2'' = map normaliseL
normaliseL :: [Double] -> [Double]                                                                        
normaliseL xs = normalise' xs (sum xs)
normalise' xs tot = map (/tot) xs

posteriorTipsLkl :: Int -> PatternAlignment -> DNode -> [[(Int,String,DNode,[Double])]]
posteriorTipsLkl numModels aln root@(DTree _ _ _ _ _ _ _) = map (map (\(i,j,k)->(i,j,k,likelihoodsPerSite aln k)))$ posteriorTipsList numModels root

posteriorTipsList :: Int -> DNode -> [[(Int,String,DNode)]]
posteriorTipsList numModels (DLeaf name dist seq partial model pl) = [map getLeaf [0..(numModels-1)]] where
                                                                        getLeaf i = (i,name,(DLeaf name dist seq partial' model (calcLeafPL partial' dist model))) where
                                                                                partial' = blankMat partial numModels i

posteriorTipsList numModels (DINode l r bl model pl) = (map (map getINodeR) $ posteriorTipsList numModels r)++ (map (map getINodeL) $ posteriorTipsList numModels l) where
                                                           getINodeL (i,name,node) = (i,name,DINode node r bl model (calcPL node r bl model))
                                                           getINodeR (i,name,node) = (i,name,DINode l node bl model (calcPL l node bl model))

posteriorTipsList numModels (DTree l m r pLs model priors pis) = (map (map getINodeR) $ posteriorTipsList numModels r) ++ (map (map getINodeM) $ posteriorTipsList numModels m) ++ (map (map getINodeL) $ posteriorTipsList numModels l) where
                                                           getINodeL (i,name,node) = (i,name,DTree node m r (calcRootPL node m r) model priors pis)
                                                           getINodeM (i,name,node) = (i,name,DTree l node r (calcRootPL l node r) model priors pis)
                                                           getINodeR (i,name,node) = (i,name,DTree l m node (calcRootPL l m node) model priors pis)
                                                                
                                                           


blankMat :: Matrix Double -> Int -> Int -> Matrix Double
blankMat matrix numClass except = fromRows filtRows where
                                        blocksize = (rows matrix) `div` numClass
                                        lower = except * blocksize
                                        upper = lower + blocksize
                                        myRows = toRows matrix
                                        zeroVect = mapVector (\x->0.0) (head myRows)
                                        filtRows = map (\(i,j) -> if (j >= lower && j < upper) then i else zeroVect) (zip myRows [0..])
                                        

likelihoods :: DNode -> [Double]
likelihoods t@(DTree _ _ _ pLs _ priors pis) = map summarise likelihoodss where
                                                summarise lkl = sum $ zipWith (*) lkl priors
                                                likelihoodss = transpose $ rawlikelihoods t

-- likelihoods per pattern, indexed by model then site
rawlikelihoods :: DNode -> [[Double]] 
rawlikelihoods (DTree _ _ _ pLs _ priors pis) = likelihoodss where
                                                likelihoodss = map toList $ map (\(pi,pL) -> pi <> pL) $ zip pis pLs
                                        

likelihoodsPerSite aln tree = likelihoodss' where
                                likelihoodss = likelihoods tree
                                mapping = mapBack aln
                                likelihoodss' = map (likelihoodss !!) mapping
                                                                

logLikelihood = logLikelihoodModel

data PatternAlignment = PatternAlignment {names :: [String], seqs::[String], columns::[String], patterns::[String], counts::[Int]}

pAlignment (ListAlignment names seqs columns) = PatternAlignment names seqs columns patterns counts where
                                                   --(patterns,counts) = unzip $ map (\x -> (head x,length x)) $ group $ sort columns
                                                   --lets keep the order (initial slow impl)
                                                   patterns = nub columns
                                                   counts = map (\x->length $ findIndices (==x) columns) patterns 
mapBack :: PatternAlignment -> [Int]
mapBack (PatternAlignment _ _ columns patterns _) = map (\x->fromJust $ findIndex (==x) patterns) columns

calcPL :: DNode -> DNode -> [Double] -> [BranchModel] -> [Matrix Double]
calcPL left right dist model = (allPLC (map (\x-> fst $ x dist) model) $ allCombine (getPL left) (getPL right)) 
calcLeafPL :: Matrix Double -> [Double] -> [BranchModel] -> [Matrix Double]
calcLeafPL tips dist model = map (\pT -> rawPartialLikelihoodCalc pT tips) $ (map (\x -> fst $ x dist) model)
calcRootPL left middle right = allCombine (getPL middle) $ allCombine (getPL left) (getPL right)

{-- gives [[[Matrix Double]]] where [branch [proc [left/right Matrix]]] --}
{-- need to map from [branch [l/r [proc Matrix]]] --}
toPBEList ((f,s,_,_,_):xs) = (map (\(i,j) -> [i,j]) $ zip f s):(toPBEList xs)
toPBEList [] = []

toPBELengths ((_,_,dist,_,_):xs) = dist:(toPBELengths xs)
toPBELengths [] = []

toPBESplits = map (\(_,_,_,a,_) -> a)
toPBEQ :: [([Matrix Double],[Matrix Double],Double,([String],[String]),[Matrix Double])] -> [[Matrix Double]]
toPBEQ = map (\(_,_,_,_,a) -> a)

--getQ node = concat $ getQ' node [] 
getQ = getQ'
getQ' (DLeaf _ bl  _  _ model _)   = (map (\x-> snd $ x bl) model) 
getQ' (DINode l r bl model _) = map (\x-> snd $ x bl) model


getPartialBranchEnds :: DNode -> [([Matrix Double],[Matrix Double],Double,([String],[String]),[Matrix Double])]
getPartialBranchEnds tree = getPartialBranchEnds' $ allRootings tree
getPartialBranchEnds' (x:xs) = (getPartialBranchEndsLeaves'' x) ++ (getPartialBranchEndsX' xs) -- first tree is special case, only do leaves
getPartialBranchEndsX' (x:xs) = (getPartialBranchEnds'' x) ++ (getPartialBranchEndsX' xs)
getPartialBranchEndsX' [] = []
getPartialBranchEnds'' (DTree left middle right pLS _ _ _) = (branchPL right,noright,head $ getBL right,splits,getQ right): remainder where
                                                                noright = allCombine (getPL left) (getPL middle) 
                                                                remainder = (calc left middle) ++ (calc middle left)
                                                                calc node other = case (node,other) of 
                                                                                       (l@(DLeaf _ bl _ leafPL _ _),o) -> [(branchPL l,allCombine (getPL o) (getPL right),head $ getBL l,((descendents o) ++ (descendents right), descendents l),getQ l)]
                                                                                       (_,o) -> []
                                                                splits = ((descendents left) ++ (descendents middle), (descendents right))

getPartialBranchEndsLeaves'' (DTree left middle right pLS _ _ _) = remainder where
                                                                noright = allCombine (getPL left) (getPL middle) 
                                                                remainder = (calc right middle left) ++ (calc left middle right) ++ (calc middle left right)
                                                                calc node other other'= case (node,other,other') of 
                                                                                       (l@(DLeaf _ bl _ leafPL _ _),o,o') -> [(branchPL l,allCombine (getPL o) (getPL o'),head $ getBL l,((descendents o) ++ (descendents o'), descendents l),getQ l)]
                                                                                       (_,_,_) -> []


branchPL :: DNode -> [Matrix Double]
branchPL (DLeaf _ _ _ leafPL _ endPL) = map (\x->leafPL) endPL
branchPL (DINode l r _ _ endPL) = map (\(x,y)-> combinePartial x y) (zip (getPL l) (getPL r))

restructData node model priors pi = restructDataMapped node (\x -> model) priors pi

restructDataMapped :: DNode -> (DNode -> [BranchModel]) -> [Double] -> [Vector Double] -> DNode 
restructDataMapped leaf@(DLeaf name dist sequence partial _ _) model priors pi  = DLeaf name dist sequence partial myModel partial' where
                                                                               partial' = calcLeafPL partial dist myModel
                                                                               myModel = model leaf

restructDataMapped inode@(DINode left right dist _ _ ) model priors pi= DINode newleft newright dist myModel partial where
                                                            partial = calcPL newleft newright dist myModel
                                                            myModel = model inode
                                                            newleft = restructDataMapped left model priors pi
                                                            newright = restructDataMapped right model priors pi

restructDataMapped tree@(DTree left middle right _ pc _ _ ) model priors pi = DTree newleft newmiddle newright partial pc priors pi where 
                                            partial = calcRootPL newleft newmiddle newright
                                            newleft = restructDataMapped left model priors pi
                                            newright = restructDataMapped right model priors pi
                                            newmiddle = restructDataMapped middle model priors pi  

removeModel (DTree l m r _ pc _ _ )  = NTree (removeModel l) (removeModel m) (removeModel r) pc
removeModel (DINode l r d _ _ ) = NINode (removeModel l) (removeModel r) d
removeModel (DLeaf name dist seq partial _ _) = NLeaf name dist seq partial



structDataN :: Int -> SeqDataType -> PatternAlignment -> Node -> NNode

structDataN hiddenClasses seqDataType pAln node = structDataN' hiddenClasses seqDataType pAln (transpose $ patterns pAln) node 
structDataN' hiddenClasses seqDataType (PatternAlignment names seqs columns patterns counts) transPat (Leaf name dist) = ans where
                                                                                                                          partial = getPartial hiddenClasses seqDataType sequence 
                                                                                                                          sequence = snd $ fromJust $ find (\(n,_) -> n==name) $ zip names transPat 
                                                                                                                          ans = NLeaf name [dist] sequence partial 


structDataN' hiddenClasses seqDataType pA transPat (INode c1 c2 dist) = NINode left right [dist] where
                                                                              left = structDataN' hiddenClasses seqDataType pA transPat c1
                                                                              right = structDataN' hiddenClasses seqDataType pA transPat c2

structDataN' hiddenClasses seqDataType pA transPat (Tree (INode l r dist) c2 ) = NTree left middle right $ counts pA where
                                                                                  left = structDataN' hiddenClasses seqDataType pA transPat l
                                                                                  right = structDataN' hiddenClasses seqDataType pA transPat r
                                                                                  middle = structDataN' hiddenClasses seqDataType pA transPat $ addDist dist c2
                                                                                  addDist d (Leaf name dist) = Leaf name (dist+d)
                                                                                  addDist d (INode l r dist) = INode l r (dist+d)

structDataN' hiddenClasses seqDataType pA transPat (Tree c2 (INode l r dist)) = structDataN' hiddenClasses seqDataType pA transPat (Tree (INode l r dist) c2)

structDataN' hiddenClasses seqDataType pA transPat (Tree (Leaf name dist) (Leaf name2 dist2))  = error "Can't handle tree with only 2 leaves"


structData hiddenClasses seqDataType pAln model transPat node = addModel (structDataN hiddenClasses seqDataType pAln node) model 



addModelNNode :: NNode -> [BranchModel] -> [Double] -> [Vector Double] -> DNode
addModelNNode (NLeaf name dist sequence partial) model _ _  = DLeaf name dist sequence partial model partial' where
                                                                                              partial' = calcLeafPL partial dist model

addModelNNode (NINode c1 c2 dist) model priors pi =  DINode left right dist model myPL where
                                                  left = addModelNNode c1 model priors pi
                                                  right = addModelNNode c2 model priors pi
                                                  myPL = calcPL left right dist model 
                                                                            
addModelNNode (NTree c1 c2 c3 pat) model priors pi = DTree left middle right myPL pat priors pi where
                                  left = addModelNNode c1 model priors pi
                                  middle = addModelNNode c2 model priors pi
                                  right= addModelNNode c3 model priors pi
                                  myPL = calcRootPL left middle right

data SeqDataType = AminoAcid | Nucleotide deriving Show
-- produce the partial likelihood at the tips
-- getPartial numHiddenStates dataType Sequence
getPartial :: Int -> SeqDataType -> String -> Matrix Double
getPartial a b c = trans $ fromRows $ getPartial' a b c
getPartial' :: Int -> SeqDataType -> String -> [Vector Double]
getPartial' _ _ [] = []
getPartial' classes AminoAcid (x:xs) = (aaPartial classes x):(getPartial' classes AminoAcid xs)
getPartial' classes Nucleotide (x:xs) = (nucPartial classes x):(getPartial' classes Nucleotide xs)

numberedHash list = foldr (\(item,i) init -> HM.insert item i init) HM.empty $ zip list [0..]

myOrderVec Nucleotide = nucOrderVec
myOrderVec AminoAcid = aaOrderVec

nucOrder = ['T','C','A','G'] --as PAML
nucOrderVec = DVec.fromList nucOrder
nucOrderHash = numberedHash nucOrder

nucPartial classes x | isGapChar x = buildVector (4*classes) (\i->1.0)
                     | otherwise = case HM.lookup x nucOrderHash of 
                                        Just index -> buildVector (4*classes) (\i -> if i`mod`4==index then 1.0 else 0.0)
                                        Nothing -> error $ "character " ++ [x] ++ " is not a nucleotide or gap"



aaOrder = ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
aaOrderVec = DVec.fromList aaOrder
aaOrderHash = numberedHash aaOrder
aaPartial classes x | isGapChar x = buildVector (20*classes) (\i->1.0) 
                    | otherwise = case HM.lookup x aaOrderHash of 
                                       Just index -> buildVector (20*classes) (\i -> if i`mod`20==index then 1.0 else 0.0) where
                                       Nothing -> error $ "character " ++ [x] ++ " is not an amino acid or gap"

                                 


quickGamma' numCat tree priors pi s aln alpha = logLikelihood dataTree  where
                                                dataTree = structData 1 AminoAcid pAln branchModel transpats tree priors (repeat pi)
                                                branchModel = (map (\r->(\x->(scaledpT r eigenS x,mat))) rates) 
                                                rates = gamma numCat alpha
                                                transpats = transpose $ patterns pAln
                                                (eigenS,mat) = quickEigen' pi s
                                                pAln = pAlignment aln
                                                patcounts = counts pAln
                                                
                                                        
quickEigen pi s = eigQ (normQ (makeQ s pi) pi) pi
quickEigen' pi s = (eigQ mat pi,mat) where
                        mat = normQ (makeQ s pi) pi

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


type BranchModel = [Double] -> (Matrix Double,Matrix Double) -- first matrix exponentiated for BL, second not

type OptModel = DNode -> [Double] -> Double
data DNode = DLeaf {dName :: String,dDistance :: [Double],sequence::String,tipLikelihoods::Matrix Double,p::[BranchModel],partialLikelihoods::[Matrix Double]} |
        DINode DNode DNode [Double] [BranchModel] [Matrix Double] | 
        DTree DNode DNode DNode [Matrix Double] [Int] [Double] [(Vector Double)]

data NNode = NLeaf String [Double] String (Matrix Double) | NINode NNode NNode [Double] | NTree NNode NNode NNode [Int]
data DataModel = DataModel {dTree::DNode, patterncounts::[Int], priors::[Double], pis::[Vector Double]}

getQMat (DTree l _ _ _ _ _ _) = getQMat l
getQMat (DLeaf _ bl _ _ m _) = map (\x-> (snd $ x bl)) m
getQMat (DINode _ _ bl m _) = map (\x-> (snd $ x bl)) m
getPi (DTree _ _ _ _ _ _ pi) = pi

getLeaves node = getLeaves' node []
getLeaves' :: DNode -> [DNode] -> [DNode]
getLeaves' (DTree l m r _ _ _ _) xs = getLeaves' r $ getLeaves' m $ getLeaves' l xs
getLeaves' (DINode l r _ _ _) xs = getLeaves' l $ getLeaves' r xs
getLeaves' leaf xs = (leaf:xs)

instance Show DNode where
        show (DTree l m r _ _ _ _) = "("++(show l)++","++(show m)++","++(show r)++");"
        show (DINode l r [d] _ _) = "(" ++ (show l)++","++(show r)++"):"++(show d)
        show (DINode l r d _ _) = "(" ++ (show l)++","++(show r)++"):"++(show d)
        show (DLeaf name [d] _ _ _ _) = name++":"++(show d)
        show (DLeaf name d _ _ _ _) = name++":"++(show d)

instance NFData DNode where
        rnf (DTree l m r mats ints doubles vectors) = rnf l `seq` rnf m `seq` rnf r `seq` rnf (map toLists mats) `seq` rnf ints `seq` rnf doubles `seq` rnf (map toList vectors) `seq` ()
        rnf (DINode l r doubles models mats) = rnf l `seq` rnf r `seq` rnf doubles `seq` rnf models `seq` rnf (map toLists mats) `seq` ()
        rnf (DLeaf string doubles s2 mat bms mats) = rnf string `seq` rnf doubles `seq` rnf s2 `seq` rnf (toLists mat) `seq` rnf bms `seq` rnf (map toLists mats) `seq` ()

instance NFData NNode where
        rnf (NLeaf s d s2 mat) = rnf s `seq` rnf d `seq` rnf s2 `seq` rnf (toLists mat)
        rnf (NINode l r bl) = rnf l `seq` rnf r `seq` rnf bl
        rnf (NTree l m r i) = rnf l `seq` rnf m `seq` rnf r `seq` rnf i

instance NFData PatternAlignment where
        rnf (PatternAlignment n s c pats counts) = rnf n `seq` rnf s `seq` rnf pats `seq` rnf c `seq` rnf counts

class AddTree a where
        addModel :: (AddTree a) => a -> [BranchModel] -> [Double] -> [Vector Double] -> DNode
        addModelFx :: (AddTree a) => a -> ([BranchModel],[Vector Double]) -> [Double] -> DNode
        addModelFx t (bm,pi) prior = addModel t bm prior pi
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
getPL (DTree _ _ _ lkl _ _ _) = lkl

getDist (DLeaf _ dist _ _ _ _) = dist
getDist (DINode _ _ dist _ _) = dist
getDist root = []

getModel (DLeaf _ _ _ _ model _) = model
getModel (DINode _ _ _ model _) = model
getModel root = error "No model for root"


goLeft (DTree (DINode l r dist model pL) middle right _ pC priors pi)  = DTree l r r' rootpL pC priors pi where
        rootpL = calcRootPL l r r'
        r' = DINode middle right dist model pL' 
        pL' = calcPL middle right dist model

goRight (DTree left middle (DINode l r dist model pL) _ pC priors pi) = DTree l' l r rootpL pC priors pi where
        rootpL = calcRootPL l' l r
        l' = DINode left middle dist model pL'
        pL' = calcPL left middle dist model

canGoLeft (DTree (DINode _ _ _ _ _) _ _ _ _ _ _ ) = True
canGoLeft _ = False
canGoMid (DTree _ (DINode _ _ _ _ _) _ _ _ _ _ ) = True
canGoMid _ = False

optBLD = optBLDx optLeftBL
optBLD0 = optBLDx (optLeftBLx 0)

optBLDFull0 :: DNode -> DNode
optBLDFull0 t = optBLDFull  (optLeftBLx 0) t
optBLDFull f tree = optBLDFull' (logLikelihood tree) 1E-05 f tree
optBLDFull' old cutoff f tree = case (new - old) of 
        x | x > cutoff -> optBLDFull' new cutoff f newtree
        otherwise      -> newtree 
        where newtree = optBLDx f tree
              new = logLikelihood newtree
        
                                
optBLDx f tree = t3 where
        t1 = swapLeftMiddle $ f $ if (canGoLeft tree) then goRight $ optBL' f (goLeft tree) else tree
        t2 = swapLeftRight $ f $ swapLeftMiddle $ f $ if (canGoLeft t1) then goRight $ optBL' f (goLeft t1) else t1
        t3 = swapLeftRight $ f $ if (canGoLeft t2) then goRight $ optBL' f (goLeft t2) else t2

optBL' f tree = ans where
        leftOptd = if (canGoLeft tree) then goRight $ optBL' f (goLeft tree) else tree
        ansL = swapLeftMiddle $ f leftOptd
        midOptd = if (canGoLeft ansL) then goRight $ optBL' f (goLeft ansL) else ansL
        ans = swapLeftMiddle $ f midOptd

allRootings tree = ans where
        midTree = swapLeftMiddle tree
        rightTree = swapLeftRight tree
        leftRootings = if (canGoLeft tree) then (tree' : allLM (tree')) else [] 
        midRootings = if (canGoLeft midTree) then (midTree' : allLM (midTree')) else [] 
        rightRootings = if (canGoLeft rightTree) then (rightTree' : allLM (rightTree')) else [] 
        tree' = goLeft tree
        midTree' = goLeft midTree
        rightTree' = goLeft rightTree
        ans  = tree : (leftRootings ++ midRootings ++ rightRootings)

getAllF2 tree = tree : ((map swapLeftMiddle (allRootings tree))  ++ (allRootings tree))

-- | return all rootings required such that each branch in the tree is given as the  
-- | the left branch in (DTree l m r _ _ _ _ ) in the *same order* as used by
-- | setBL and getBL
-- | setLeftBL can then be used to set the branch lengths
getAllF :: DNode -> [DNode]
getAllF tree = ans where
        t1 = if (canGoLeft tree) then getAllF' (goLeft tree) else [] --leftRootings
        t1' = swapLeftMiddle tree
        t2 = if (canGoLeft t1') then getAllF' (goLeft t1') else [] -- midRootings 
        t2' = swapLeftRight tree
        t3 = if (canGoLeft t2') then getAllF' (goLeft t2') else [] -- rightRootings
        ans = concat [(tree:t1),(t1':t2),(t2':t3)]
        

getAllF' tree = (concat [(tree:t1),(t1':t2)]) where
        t1 = if (canGoLeft tree) then getAllF' (goLeft tree) else [] -- leftRootings
        t1' = swapLeftMiddle tree 
        t2 = if (canGoLeft t1') then getAllF' (goLeft t1') else [] -- midRootings



getLeftSplit node = getLeftSplit' node []
getLeftSplit' (DLeaf name bl _ _ _ _) spls = [name]:spls
getLeftSplit' (DINode l r bl _ _) spls = ((descendents l) ++ (descendents r)):(getLeftSplit' l (getLeftSplit' r spls))
getLeftSplit' (DTree l m r _ _ _ _ ) spls = getLeftSplit' l ( getLeftSplit' m ( getLeftSplit' r spls))

getBL node = concat $ getBL' node []

getPriors (DTree _ _ _ _ _ priors _ ) = priors

getBL' (DLeaf _ bl _ _ _ _) bls = bl : bls
getBL' (DINode l r bl _ _) bls = bl : (getBL' l (getBL' r bls))
getBL' (DTree l m r _ _ _ _ ) bls = (getBL' l (getBL' m (getBL' r bls)))

allLM tree = a ++ b where
               a = allLeftRootings tree
               b = allMidRootings tree

allMidRootings tree | canGoMid tree = t2 : (allLM t2) where
                                      t2 = goLeft $ swapLeftMiddle tree
allMidRootings tree = []

allLeftRootings tree | canGoLeft tree = t2 : (allLM t2) where
                                        t2 = goLeft tree
allLeftRootings tree = []

descendents (DTree l m r _ _ _ _) = (descendents l) ++ (descendents m ) ++ (descendents r)
descendents (DLeaf name _ _ _ _ _) = [name] 
descendents (DINode l r _ _ _) = (descendents l) ++ (descendents r)

nonRootNodes (DTree l m r _ _ _ _) = (nonRootNodes l) ++ (nonRootNodes m) ++ (nonRootNodes r)
nonRootNodes node@(DINode l r _ _ _) = node : ((nonRootNodes l) ++ (nonRootNodes r))
nonRootNodes node@(DLeaf _ _ _ _ _ _) = [node]

--setLeftBL (DTree l m r _ _ _ _ ) x | trace ((show x ) ++ " (" ++ (descendents l) ++ " : " ++ (descendents m) ++ " : " ++ (descendents r) ++ ")") False = undefined

setLeftBL0 node x = setLeftBL node x' where
        x' = x:(tail (getLeftBL node))

setLeftBL1 node x = setLeftBL node x' where
        x' = (head (getLeftBL node)):x'

setLeftBL (DTree (DLeaf name dist seq tip model _) m r _ pC priors pi) x = (DTree newleft m r pl' pC priors pi) where
        partial' = calcLeafPL tip x model
        newleft = DLeaf name x seq tip model partial'
        pl' = calcRootPL newleft m r 


setLeftBL (DTree (DINode l1 r1 dist model _) m r _ pC priors pi) x = (DTree newleft m r pl' pC priors pi) where
        partial' = calcPL l1 r1 x model
        newleft = DINode l1 r1 x model partial'
        pl' = calcRootPL newleft m r 

--optLeftBL tree  = traceShow ("OptBL " ++ (show best) ++ "->" ++ (show $ logLikelihood $ setLeftBL tree best)) $ setLeftBL tree best where
optLeftBL tree  = setLeftBL tree best where
                        (DTree l _ _ _ _ _ _ ) = tree
                        best = case l of 
                                (DINode _ _ [param] _ _) -> getLeftBL $ optLeftBLx 0 tree
                                (DLeaf _ [param] _ _ _ _) -> getLeftBL $ optLeftBLx 0 tree
                                (DLeaf _ params _ _ _ _) -> fst $ maximize NMSimplex2 1E-4 1000 (map (\i->0.05) params) (boundedBLfunc (\x -> logLikelihood $ setLeftBL tree x)) params    
                                (DINode _ _ params _ _) -> fst $ maximize NMSimplex2 1E-4 1000 (map (\i->0.05) params) (boundedBLfunc (\x -> logLikelihood $ setLeftBL tree x)) params    

getLeftBL (DTree (DINode _ _ x _ _) _ _ _ _ _ _) = x
getLeftBL (DTree (DLeaf _ x _ _ _ _) _ _ _ _ _ _) = x
#ifdef Debug
sensibleBrent tol ftol f hardL hardU startB startF = mytrace  ("BRENT " ++ (show boundL') ++ " " ++ (show boundU') ++ " " ++ (show startB')) $ brent tol ftol f boundL boundU startB' fBoundL fBoundU startF' where
#else
sensibleBrent tol ftol f hardL hardU startB startF = brent tol ftol f boundL boundU startB' fBoundL fBoundU startF' where
#endif
        boundL' = max hardL (startB/2)
        boundU' = case min hardU (startB*2) of
                      x | x <= boundL' -> hardU 
                        | otherwise    -> x 
        (startB',startF') = case startB of
                                x | x > boundL' && x < boundU' -> (startB,startF)
                                  | otherwise                  -> ((boundU' + boundL')/2,f $ (boundU'+boundL')/2)
        fBoundU' = f boundU'
        fBoundL' = f boundL'
        (boundL,fBoundL) = case (boundL',fBoundL') of 
                                (a,b) | b > startF -> (boundL',fBoundL')
                                      | otherwise  -> (hardL, f hardL)
        (boundU,fBoundU) = case (boundU',fBoundU') of 
                                (a,b) | b > startF -> (boundU',fBoundU')
                                      | otherwise  -> (hardU, f hardU)


--optLeftBLx val tree  = traceShow ("OptBL " ++ (show startBL) ++ " -> " ++ (show best) ++ " -> " ++ (show $ logLikelihood $ setLeftBL tree best)) $ bestTree where
optLeftBLx val tree  =  bestTree where
                                startBL = getLeftBL tree
                                startB = head $ drop val startBL
                                startLL = logLikelihood tree
                                f = (loggedFunc $ invert $ (\x -> logLikelihood $ setLeftBL tree (replace val [x] startBL)))
                                b= sensibleBrent 1E-3 1E-6 f 1E-6 50.0 startB (f startB)
                                best = replace val [b] startBL
                                bestTree' = setLeftBL tree best 
                                --protect from golden section screw-ups
                                bestTree = case (logLikelihood bestTree') of 
                                                x | x < startLL -> tree
                                                  | otherwise -> bestTree'

boundedBLfunc = boundedFunction (-1E20) (repeat $ Just 0.0) (repeat Nothing)

swapLeftMiddle (DTree l m r pL pC priors pi) = DTree m l r pL pC priors pi
swapLeftRight (DTree l m r pL pC priors pi) = DTree r m l pL pC priors pi


getSplits node = map invertSplit $ getSplits' node [] where
                        alldescendents = sort $ descendents node
                        invertSplit s = ((sort s), filter (\x -> not $ x `elem` s) alldescendents) 

getSplits' (DTree l m r _ _ _ _ ) splits = (getSplits' l (getSplits' m (getSplits' r splits)))
getSplits' (DINode l r _ _ _ ) splits = (descendents l ++ descendents r):(getSplits' l (getSplits' r splits))
getSplits' (DLeaf name _ _ _ _ _) splits = [name]:splits




-- | Converts a mapping ('DNode'->'Int') into ['Int'], in the correct order to iterate over the branches
-- | using setBL and variants
getLinearMap mapping (DTree l m r _ _ _ _ ) ints = (getLinearMap mapping l (getLinearMap mapping m (getLinearMap mapping r ints)))
getLinearMap mapping node@(DINode l r _ _ _) ints = (mapping node) : (getLinearMap mapping l (getLinearMap mapping r ints))
getLinearMap mapping leaf ints = (mapping leaf) : ints


{-setBL :: [Double] -> DNode -> DNode -}
{-setBL bls node = fst $ setBL' reformattedBLs node where-}
        {-oldBLs = getBL' node []-}
        {-reformattedBLs = reformat bls oldBLs-}
        {-reformat blList (x:xs) = (taken : reformat remainder xs) where-}
                {-(taken,remainder) = splitAt (length x) blList-}

setBL' x = setBLX' 0 x

setLeftBLMapped i (DTree l m r _ pC priors pi) mapping = DTree left m r partial pC priors pi where
        left = setBLMapped i l mapping
        partial = calcRootPL left m r

setBLMapped i (DTree l m r _ pC priors pi) mapping = DTree left middle right partial pC priors pi where
        left = setBLMapped i l mapping
        middle = setBLMapped i m mapping
        right = setBLMapped i r mapping
        partial = calcRootPL left middle right

setBLMapped i (DINode l r blstart mats pl) mapping = DINode left right newbl mats partial where
        left = setBLMapped i l mapping
        right = setBLMapped i r mapping
        partial = calcPL left right (newbl) mats
        initBL = take i blstart
        endBL = mapping (DINode left right blstart mats pl)
        newbl = initBL ++ endBL

setBLMapped i (DLeaf a blstart b tips model pl) mapping = DLeaf a newbl b tips model partial where
        partial = calcLeafPL tips newbl model
        initBL = take i blstart
        endBL = mapping (DLeaf a blstart b tips model pl)
        newbl = initBL ++ endBL

leftSplit (DTree l m r _ _ _ _ ) = (descendents l, (descendents m) ++ (descendents r))


setBLX' i bls (DTree l m r pl pC priors pi) = ((DTree left middle right partial pC priors pi), remainder3) where
                        (left,remainder) = setBLX' i bls l
                        (middle,remainder2) = setBLX' i remainder m
                        (right,remainder3) = setBLX' i remainder2 r
                        partial = calcRootPL left middle right

setBLX' i (bl2:bls) (DINode l r blstart mats pl) = ((DINode left right newbl mats partial), remainder2) where
                             (left,remainder) = setBLX' i bls l
                             (right,remainder2) = setBLX' i remainder r
                             partial = calcPL left right (newbl) mats
                             newbl = replace i bl2 blstart

setBLX' i (bl2:bls) (DLeaf a (blstart) b tips model _) = ((DLeaf a newbl b tips model pl),bls) where
                                      pl = calcLeafPL tips newbl model 
                                      newbl = replace i bl2 blstart

setBLX' _ [] _ = error "Insuffiencent branch lengths specified"


replace :: Int -> [a] -> [a] -> [a]
replace 0 [x] (s:ss) = (x:ss)
replace i x start =  (take (i) start) ++ x ++ (drop (i + (length x)) start)

listpartition (c:counts) start ans = listpartition counts b (a:ans) where
                                 (a,b) = splitAt c start
listpartition [] [] ans = reverse ans
listpartition [] start ans = listpartition [] [] (start:ans)


setBLflat i node = setBLflat' i node $ map (\x -> (length x) - i) (getBL' node [])
setBLflat' i node counts bls = setBLX' i (listpartition counts bls []) node
                        

zeroQMat size = diag $ constant 0.0 size

qdLkl numCat priors dataType aln tree modelF params = logLikelihood $ addModelFx (structDataN numCat dataType (pAlignment aln) tree) (modelF params) priors

type ModelF = [Double] -> ([BranchModel],[Vector Double])
basicModel :: Matrix Double -> Vector Double -> ModelF
basicModel s pi _ = ([model],[pi]) where
                    model = (\x -> (standardpT e x,m))
                    (e,m) = quickEigen' pi s

quickLkl aln tree pi s = qdLkl 1 [1.0] AminoAcid aln tree (basicModel s pi) []
quickGamma numCat alpha aln tree pi s = qdLkl 1 (flatPriors numCat) AminoAcid aln tree (gammaModel numCat s pi) [alpha]

gammaModel numCat (sF,sN) (piF,piN) (alpha:xs) = (models,replicate numCat pi) where
                                models = map (\r -> (\x ->(scaledpT r eigenS x,mat))) $ gamma numCat alpha
                                (eigenS,mat) = quickEigen' pi s
                                s = sF (take sN xs)
                                pi = piF (take piN (drop sN xs))

gammaModelQ numCat (sF,sN) (piF,piN) (alpha:xs) = map (\r -> setRate r initMat pi) $ gamma numCat alpha where
                                initMat = makeQ s pi
                                s = sF (take sN xs)
                                pi = piF (take piN (drop sN xs))
                                


thmmModel :: Int -> ([Double]->Matrix Double,Int) -> ([Double] -> Vector Double,Int) -> [Double] -> ([BranchModel],[Vector Double])
thmmModel numCat (sF,sN) (piF,piN) (priorZero:alpha:sigma:xs) = ([thmm numCat pi s priors alpha sigma],[fullPi]) where
                               priors =  (replicate (numCat -1) ((1.0-priorZero)/(fromIntegral (numCat-1))) ) ++ [priorZero]
                               s = sF (take sN xs)
                               pi = piF (take piN (drop sN xs))
                               fullPi = makeFullPi priors pi

{-- TODO refactor this out --}
thmmModelQ numCat (sF,sN) (piF,piN) (priorZero:alpha:sigma:xs) = (thmmQ numCat pi s priors alpha sigma) where
                                   priors = (replicate (numCat -1) ((1.0-priorZero)/(fromIntegral (numCat-1))) ) ++ [priorZero]
                                   fullPi = makeFullPi priors pi
                                   s = sF (take sN xs)
                                   pi = piF (take piN (drop sN xs))


thmmPerBranchModel numCat s pi [priorZero,alpha] = ([thmmPerBranch numCat pi s priors alpha],[fullPi]) where 
                        priors = (replicate (numCat -1) ((1.0-priorZero)/(fromIntegral (numCat-1))) ) ++ [priorZero]
                        fullPi = makeFullPi priors pi
thmmPerBranchModel numCat s pi list = error $ "Fail " ++ (show list)

quickThmm numCat aln tree pi s [priorZero,alpha,sigma] = qdLkl numCat [1.0] AminoAcid aln tree (thmmModel numCat s pi) [priorZero,alpha,sigma] 
                                                        
thmm :: Int -> Vector Double -> Matrix Double -> [Double] -> Double -> Double -> BranchModel
thmm numCat pi s priors alpha sigma = (\x->(standardpT eigM x,mat)) where
                                               eigM = eigQ mat fullPi
                                               mat = thmmQ numCat pi s priors alpha sigma
                                               fullPi = makeFullPi priors pi


     
--thmmQ :: Int -> Vector Double -> Matrix Double -> [Double] -> Double -> Double -> (Matrix Double,(Matrix Double, Vector Double, Matrix Double))
thmmQ :: Int -> Vector Double -> Matrix Double -> [Double] -> Double -> Double -> Matrix Double
thmmQ numCat pi s priors alpha sigma = (matDr `kronecker` matM) `add` (matG `kronecker` matIm) where 
                                        matM = normQ (makeQ s pi) pi
                                        matDr = diag $ fromList $ map (*factor) ((gamma (numCat -1) alpha) ++ [0.0])
                                        factor = 1.0 / (1.0 - (last priors))                                                                                                                                                                  
                                        matG = makeQ sigmaMat $ fromList priors
                                        sigmaMat = diagRect sigma (fromList $ replicate numCat 0.0) numCat numCat 
                                        matIm = diag pi 

thmmPerBranch :: Int -> Vector Double -> Matrix Double -> [Double] -> Double -> [Double] -> (Matrix Double,Matrix Double)
thmmPerBranch numCat pi s priors alpha [branchLength,sigma] = thmm numCat pi s priors alpha sigma [branchLength]
thmmPerBranch numCat pi s priors alpha paramList | mytrace ("thmmpb " ++ (show paramList)) False = undefined


loggedFuncGeneric :: (Show a) => (a->b) -> (a->b)
loggedFuncGeneric f x = trace (show x) (f x)

loggedFunc :: (Show a) => (a -> Double) -> (a -> Double)
loggedFunc f | mytrace "LOGGING ON" True = f2 where
               f2 x = ans3 where
                      ans3 | mytrace (show x) True = ans2
                      ans2 | mytrace ((show x) ++ " -> " ++ (show ans)) True = ans
                      ans = f x
                    
maximize method pre maxiter size f = minimize method pre maxiter size (invert f) 

stepSize = 1E-4

optBL model t' params priors cutoff = mytraceShow optTree optTree where
        tree  = addModelFx t' (model params) priors
        optTree = optBLD $ optBLD tree

optBLx :: (AddTree t) => Int ->  ModelF -> t -> [Double] -> [Double] -> Double -> DNode
optBLx val model t' params priors cutoff = mytraceShow optTree optTree where
                tree  = addModelFx t' (model params) priors
                optTree = optBLD $ optBLD tree
 

getSingleDimensional :: ([Double]->a) -> [Double] -> [Double->a]
getSingleDimensional f x = map (\i x2 -> f (replace i [x2] x)) $ map fst $ zip [0..] x

--getFunc1 priors model tree params = logLikelihood $ getFuncT1 priors model tree params

-- | function that sets only parameters across a tree
getFuncT1 :: [Double] -> ModelF -> DNode -> [Double] -> (DNode,[Double->DNode])
getFuncT1 priors model tree params = (newtree,funcs) where
        funcs = getSingleDimensional f params
        newtree = f params where
        f x = addModelFx tree (model x) priors

-- | function that sets bsParams then parameters 
-- | bsParams are set by a mapping (DNode -> Int) 
getFuncT1A  :: [Double] -> (DNode -> Int) -> (Int,Int) -> ModelF -> DNode -> [Double] -> (DNode,[Double->DNode])
getFuncT1A priors mapping (0,_) model tree params = getFuncT1 priors model tree params
getFuncT1A priors mapping numBSParam model tree params = (newTree,funcs) where
        f = getFuncT1A' priors mapping numBSParam model tree 
        newTree = f params
        funcs = getSingleDimensional (loggedFuncGeneric f) params

getFuncT1A' priors mapping (paramsPerBranch,paramCats) model tree params = newtree where
        modelx = model params'
        newtree = addModelFx t2 modelx priors where
                t2 = setBLMapped 1 tree getParamF 
        perBranchParams = take paramCats $ splitLists paramsPerBranch params
        params' = drop (paramsPerBranch * paramCats) params
        getParamF node = (perBranchParams !! (mapping node))

--getParamsT1 params tree = params

-- branchLengths, then other parameters
getFuncT2 :: [Double] -> ModelF -> DNode -> [Double] -> (DNode,[Double->DNode])
getFuncT2 priors model tree params = (newtree,funcs) where
        funcs = funcsPerBranch ++ (getSingleDimensional f params2)
        newtree = f params2
        f x = addModelFx t2 (model x) priors
        (t2,params2') = setBLX' 0 (map (\x->[x]) params) tree
        params2 = map head params2'
        rootedTrees = getAllF t2
        funcsPerBranch = map setLeftBL0 rootedTrees

--getParamsT2 params tree = (map head $ getBL' tree []) ++  params

-- branchLengths, then bsParams, then other params
getFuncT3 priors mapping (0,_) model tree params = getFuncT2 priors model tree params
getFuncT3 priors mapping numBSParam model tree params | mytrace ("getFuncT3 " ++ (show params)) True = (newtree,funcs) where
        (tree',remainder) = setBL' (map (:[]) params) tree
        remainder' = concat remainder
        (newtree,funcsPerParam) = getFuncT1A priors mapping numBSParam model tree' remainder'
        funcs = map setLeftBL0 rootedTrees
        rootedTrees = getAllF newtree
        

--getParamsT3 = getParamsT2 


splitLists i [] = []
splitLists i x@(z:zs) = y:(splitLists i ys) where (y,ys) = splitAt i x
                 

--optBSParamsBLIO method numBSParam mapping initialStepSize = optWithBSIO' method [] 1E-2 numBSParam (Just mapping) initialStepSize (map (/100.0) initialStepSize)
optBSParamsBL cutoff method numBSParam mapping initialStepSize = optWithBSIO' method [] cutoff numBSParam (Just mapping) initialStepSize (map (/100.0) initialStepSize)

getDeriv (curr,f,x,l,u) = getDeriv' 1E-6 curr f x l u

--getDeriv' stepSize curr f x l u | trace ("getting deriv for " ++ (show x) ++ " " ++ (show [l,u])) False = undefined 
getDeriv' stepSize curr f x l u = case (u,l) of 
                                  (Just a, Just b) -> go (min a (x+stepSize)) (max b (x-stepSize))
                                  (Just a, Nothing) -> go (min a (x+stepSize)) (x-stepSize)
                                  (Nothing,Just b) -> go (x+stepSize) (max b (x-stepSize))
                                  (Nothing,Nothing) -> go (x+stepSize) (x-stepSize) 
                                  where
                                        go a b = ((f a) - (f b))/(a-b)

data IterType = Start | Stop | Full | BL | Params deriving (Eq)
instance Show IterType where
        show Start = "Start"
        show Stop = "Finish"
        show Full = "Params/Branch"
        show BL = "Branch"
        show Params = "Params"
optWithBSIO' ::  NLOptMethod -> [(Double,IterType)] -> Double -> (Int,Int) -> (Maybe (DNode -> Int)) -> [Double] -> [Double] -> [Maybe Double] -> [Maybe Double] -> [Double] -> ModelF -> DNode -> [Double] -> [(DNode,[Double],IterType)]
optWithBSIO' method iterations cutoff numBSParam mapping stepSize limitStepSize lower upper priors model tree startParams = ans where
     ans = case iterations of
                x | mytrace ("Iterations " ++ show iterations ++ (show cutoff)) False -> undefined
    --            ((x,t):(x',t'):(x'',t''):(x''',t'''):xs) | ("Checking " ++ (show (x-x''')) ++ " " ++ (show cutoff)) False -> undefined
                ((x,t):(x',t'):(x'',t''):(x''',t'''):xs) | (x-x''' < cutoff) && (cutoff > 0.05 || (t'''==Full || t''==Full || t'==Full ||  t==Full)) -> [(tree,startParams,Stop)] --stop
    --            ((x,t):(x',t'):(x'',t''):(x''',t'''):xs) | (x-x''' < cutoff) && (t'''==Full || t==Full) -> optWithBSIO' method iterations cutoff numBSParam mapping (incrementStepSize stepSize) limitStepSize lower upper priors model tree startParams
                --((x,t):(x',t'):(x'',t''):(x''',t'''):xs) | (x-x''' < cutoff) && stepSizeMet -> return (tree,startParams) --stop
                --((x,t):(x',t'):(x'',t''):(x''',t'''):xs) | (x-x''' < cutoff) -> optWithBSIO' method iterations cutoff numBSParam mapping (incrementStepSize stepSize) limitStepSize lower upper priors model tree startParams
                list@((x,BL):[]) -> (tree',bestParams',Params) : optWithBSIO' method ((lkl,Params):list) cutoff numBSParam mapping stepSize limitStepSize lower upper priors model tree' bestParams' where
                                --initial param opt, use big steps
                                (bestParams',_) = traceShow (map (*10) stepSize) $ bobyqa (map (*10) stepSize) 1E-2 (enforceBounds lower upper startParams) f lower upper 
                                tree' = fst $ getFuncT1A priors (fromJust mapping) numBSParam model tree bestParams'
                                lkl = logLikelihood tree'
                                myFunc = getFuncT1A priors (fromJust mapping) numBSParam model tree
                                f x = mytrace ((show x) ++ " -> " ++ (show ans) ++ "\n") (ans,Just derivs) 
                                        where
                                          (outTree,outFuncs) = myFunc x
                                          ans = mytrace ((show x) ++ "...") logLikelihood outTree
                                          derivs = map getDeriv $ zip5 (repeat ans) (map (logLikelihood .) outFuncs) x lower upper
                list@((x,BL):xs) -> (tree',bestParams',Params) : optWithBSIO' method ((lkl,Params):list) cutoff numBSParam mapping stepSize limitStepSize lower upper priors model tree' bestParams' where -- Opt Params
                                bestParams2 = case (tail startParams) of 
                                                        [] -> [sensibleBrent 1E-3 1E-6 f2 upper' lower' (head startParams) (f2 $ head startParams)] where
                                                                        upper' = fromMaybe (-(1E10)) (head lower)
                                                                        lower' = fromMaybe 1E10 (head upper)
                                                                        f2 x = -(fst $ f [x])
                                                        _  -> fst $ bobyqa stepSize 1E-6 (enforceBounds lower upper startParams) f lower upper 
                                tree2 = fst $ getFuncT1A priors (fromJust mapping) numBSParam model tree bestParams2
                                lkl = logLikelihood tree2
                                (tree',bestParams') = case lkl of 
                                                        n | n > x     -> (tree2,bestParams2)
                                                          | otherwise -> (tree,startParams) --no improvement!!
                                myFunc = getFuncT1A priors (fromJust mapping) numBSParam model tree
                                f x = mytrace ((show x) ++ " -> " ++ (show ans) ++ "\n") (ans,Just derivs) where
                                        (outTree,outFuncs) = myFunc x
                                        ans = logLikelihood outTree
                                        derivs = map getDeriv $ zip5 (repeat ans) (map (logLikelihood .) outFuncs) x lower upper
                list@((x,t):(x',t'):_:_:xs) | (x-x' < (cutoff*10)) && (t /= Full) -> trace "FULL" $  (tree',bestParams',Full) : optWithBSIO' method ((lkl,Full):list) cutoff numBSParam mapping stepSize limitStepSize lower upper priors model tree' (dropBL bestParams') where --full Opt
                                stepSize' = (map (*0.1) myBL)  ++ (map (/2) stepSize)
                           {-let startParams' = enforceBounds lower' upper' (addBL startParams)-}
                                startParams' = addBL startParams
                                (bestParams',_) = method stepSize' 1E-6 startParams' f lower' upper' --f startParams'
                            --let (bestParams',err) = output
                                (tree',_) = getFuncT3 priors (fromJust mapping) numBSParam model tree bestParams'
                                lkl = logLikelihood tree'
                                myFunc = getFuncT3 priors (fromJust mapping) numBSParam model tree
                                f x = mytrace ((show x) ++ " -> " ++ (show ans)) (ans,Just derivs) where
                                        (outTree,outFuncs) = myFunc x
                                        ans = logLikelihood outTree
                                        derivs = map getDeriv $ zip5 (repeat ans) (map (logLikelihood .) outFuncs) x lower' upper'
                                myBL = map head $ getBL' tree [] 
                                numBL = length myBL
                                lower' = (replicate numBL $ Just 1E-6) ++ lower
                                upper' = (replicate numBL $ Just 50.0) ++ upper
                                addBL p = myBL ++ startParams
                                dropBL p = drop numBL p
                list@[] -> (startTree,startParams,Start):(tree',params',BL) : optWithBSIO' method ((lkl,BL):list) cutoff numBSParam mapping stepSize limitStepSize lower upper priors model tree' startParams where -- Opt Params
                                startTree = fst $ getFuncT1A priors (fromJust mapping) numBSParam model tree startParams                                            
                                tree' = optBLDFull' (logLikelihood startTree) 1E-2 (optLeftBLx 0) startTree                                       
                                params' = drop ((fst numBSParam) * (snd numBSParam)) startParams
                                lkl = logLikelihood $ fst $ getFuncT1A priors (fromJust mapping) numBSParam model tree' startParams
                list -> (tree',params',BL) : optWithBSIO' method ((lkl,BL):list) cutoff numBSParam mapping stepSize limitStepSize lower upper priors model tree' startParams where -- Opt Params
                                startTree = fst $ getFuncT1A priors (fromJust mapping) numBSParam model tree startParams                                            
                                tree' = optBLDFull0 startTree                                       
                                params' = drop ((fst numBSParam) * (snd numBSParam)) startParams
                                lkl = logLikelihood $ fst $ getFuncT1A priors (fromJust mapping) numBSParam model tree' startParams
     --incrementStepSize ss = map (\(x,y) -> max x y) $ zip limitStepSize $ map (/10) ss
     

 
enforceBounds l u p = enforce (>) u $ enforce (<) l p where
        enforce f (x:xs) (y:ys) = (enforce' f x y):(enforce f xs ys)
        enforce f [] [] = []
        enforce' f x y = case x of
                        Nothing -> y
                        Just x' -> if (f y x') then x' else y


--TODO eradicate this function!
dummyTree :: (AddTree t) => t -> DNode 
dummyTree t = addModelFx t (basicModel wagS wagPi []) [1.0]

makeMapping :: (AddTree t) => (([String],[String]) -> a) -> t -> (DNode -> a) 
makeMapping splitmap t = lookupF where
        splits = getSplits $ dummyTree t
        splitList = zip splits $ map splitmap splits 
        eitherSideList = splitList >>= (\((a,b),c)->[(a,c),(b,c)]) 
        lookupF node = case find (\(x,a) -> (sort x) == (sort $ descendents node)) eitherSideList of
                            Just (_,y) -> y
                            Nothing -> error $ "Can't find in " ++ (show (sort $ descendents node)) ++ "\nin " ++ (show $ map fst eitherSideList) 


traceX x = traceShow x x
traceXP p x = trace (p ++ " " ++ (show x)) x


genList :: StdGen -> [StdGen]
genList stdGen = first:remainder where
        (first,next) = split stdGen
        remainder = genList next

draw :: Matrix Double -> Int -> Double -> Int
--draw mat row prob | trace ((show mat) ++ "\n" ++ (show row) ++ "\n" ++ (show prob)) False = undefined
draw mat row prob  =  draw' 0 mat row prob
draw' :: Int -> Matrix Double -> Int -> Double -> Int
--draw' i mat row prob | trace ((show i) ++ " " ++  (show row) ++ " " ++ (show prob)) False = undefined
draw' i mat row prob = case (prob - (mat @@> (row,i))) of 
                          x | x <= 0.0 -> i
                            | otherwise -> draw' (i+1) mat row x

drawVec vec prob = drawVec' 0 vec prob
drawVec' i vec prob = case (prob - vec @> i) of
                          x | x <= 0.0 -> i
                            | otherwise -> drawVec' (i+1) vec x

makeSimulatedAlignment' dt stdGen t 0 = []
makeSimulatedAlignment' dt stdGen t i = (column,stdGen') : (makeSimulatedAlignment' dt stdGen' t (i-1)) where
        (column,stdGen') = makeSimulatedColumnX dt t stdGen

overlayGaps from@(ListAlignment names seqs _) to@(ListAlignment names' seqs' _) = quickListAlignment sortedNames' sortedSeqs'' where
        sortedSeqs = snd $ unzip $ sortBy (comparing snd) $ zip names seqs
        (sortedNames',sortedSeqs') = unzip $ sortBy (comparing snd) $ zip names' seqs'
        sortedSeqs'' = overlayGaps' sortedSeqs sortedSeqs'
        overlayGaps' x y = map overlayGaps'' $ zip x y
        overlayGaps'' ((x:xs),(y:ys))=case x of 
                                           x | isGapChar x -> x:(overlayGaps'' (xs,ys))
                                             | otherwise -> y:(overlayGaps'' (xs,ys))
        overlayGaps'' _ = []

makeSimulatedAlignmentWithGaps dt stdGen t oldAln@(ListAlignment n s c) = overlayGaps oldAln (makeSimulatedAlignment dt stdGen t (length $ c))

makeSimulatedAlignment dt stdGen t i = quickListAlignment names seqs where
        raw = makeSimulatedAlignment' dt stdGen t i
        columns = map fst raw
        names = map snd $ head columns
        seqs = transpose $ map (map fst) columns


makeSimulatedColumn :: SeqDataType -> DNode -> StdGen -> Int -> Int -> ([(Char,String)],StdGen)
makeSimulatedColumn dt (DLeaf name dist _ _ modelList _) gen myRow modelIndex = ans where
        ans = ([(myBase,name)],gen')
        (rand,gen') = randomR (0.0,1.0) gen 
        pT = (fst $ (modelList !! modelIndex) dist)
        orderVec = myOrderVec dt
        myBase =  orderVec ! ((draw pT myRow rand) `mod` (dataSize dt))

makeSimulatedColumn dt (DINode l r dist modelList _) gen myRow modelIndex = ans where
        ans = (myBase,genR)
        (rand,gen') = randomR (0.0,1.0) gen
        pT = (fst $ (modelList !! modelIndex) dist)
        myRow' = draw pT myRow rand
        (mybasesl,genL) = makeSimulatedColumn dt l gen' myRow' modelIndex
        (mybasesr,genR) = makeSimulatedColumn dt r genL myRow' modelIndex
        myBase = mybasesl ++ mybasesr


makeSimulatedColumn dt (DTree l m r _ _ priors pis) gen myRow modelIndex = ans where
        ans = (myBase,genR)
        (mybasesl,genL) = makeSimulatedColumn dt l gen myRow modelIndex
        (mybasesm,genM) = makeSimulatedColumn dt m genL myRow modelIndex
        (mybasesr,genR) = makeSimulatedColumn dt r genM myRow modelIndex
        myBase = mybasesl ++ mybasesm ++ mybasesr
        

makeSimulatedColumnX dt tree@(DTree l m r _ _ priors pis) gen = makeSimulatedColumn dt tree gen'' myRow modelIndex where
        (rand,gen') = randomR (0.0,1.0) gen
        (rand2,gen'') = randomR (0.0,1.0) gen'
        myRow = drawVec (pis !! modelIndex) rand
        modelIndex = drawFromDist priors rand2 

drawFromDist pris = drawFromDist' $ map nonNeg pris where
        nonNeg x | x < 0.0 = 0.0
        nonNeg x = x

drawFromDist' pris | (sum pris) < 1.0 = drawFromDist' normPris where
        normPris = map (/total) pris
        total = sum pris

drawFromDist' pris = drawFromDist'' 0 pris

drawFromDist'' :: Int -> [Double] -> Double -> Int
--drawFromDist' c pris p | trace ((show c) ++ " " ++ (show pris) ++ " " ++ (show p)) False = undefined
drawFromDist'' c (pri:pris) p = if (p <= pri) then c else (drawFromDist'' (c+1) pris (p-pri))


annotateTreeWith :: [(([String],[String]),Double)] -> DNode -> DNode 
annotateTreeWith start = annotateTreeWith' mapping where
                                listMap :: [([String],Double)] = concatMap (\((a,b),c) -> [(sort a,c),(sort b,c)]) start
                                mapping x = snd $ fromJust $ find ((==(sort x)) . fst) listMap
                   

annotateTreeWith' :: ([String]->Double) -> DNode -> DNode 
annotateTreeWith' mapping (DTree l m r models patcounts priors pis) = DTree (an l) (an m) (an r) models patcounts priors pis where
                                                                                                 an = annotateTreeWith' mapping

annotateTreeWith' mapping n@(DINode l r bl models mats) = DINode (an l) (an r) (bl ++ [(mapping $ descendents n)]) models mats where
                                                                an = annotateTreeWith' mapping

annotateTreeWith' mapping n@(DLeaf name bl seq tip models partial) = DLeaf name (bl ++ [(mapping [name])]) seq tip models partial

annotateTreeWithNumberSwitchesSigma dt sigma (DTree l m r models patcounts priors pis) = DTree (an l) (an m) (an r) models patcounts priors pis where
                                                                            an = annotateTreeWithNumberSwitchesSigma' dt sigma priors pis

annotateTreeWithNumberSwitchesSigma' dt sigma priors pis (DINode l r (bl:[]) models mats) = DINode (an l) (an r) [bl,sigma,switchbl] models mats where
                                                                                        switchbl = calcSwitchBL dt priors pis models [bl,sigma]
                                                                                        an = annotateTreeWithNumberSwitchesSigma' dt sigma priors pis

annotateTreeWithNumberSwitchesSigma' dt sigma priors pis (DLeaf name (bl:[]) seq tip models partial) = DLeaf name [bl,sigma,switchbl] seq tip models partial where
                                                                                        switchbl = calcSwitchBL dt priors pis models [bl,sigma]

annotateTreeWithNumberSwitches dt (DTree l m r models patcounts priors pis) = DTree (an l) (an m) (an r) models patcounts priors pis where
                                                                            an = annotateTreeWithNumberSwitches' dt priors pis
annotateTreeWithNumberSwitches' dt priors pis (DINode l r (bl:sigma:[]) models mats) = DINode (an l) (an r) [bl,sigma,switchbl] models mats where
                                                                                        switchbl = calcSwitchBL dt priors pis models [bl,sigma]
                                                                                        an = annotateTreeWithNumberSwitches' dt priors pis

annotateTreeWithNumberSwitches' dt priors pis (DLeaf name (bl:sigma:[]) seq tip models partial) = DLeaf name [bl,sigma,switchbl] seq tip models partial where
                                                                                        switchbl = calcSwitchBL dt priors pis models [bl,sigma]


calcSwitchBL dt priors pis models (bl:sigma:[]) = (*) bl $ sum $ zipWith (*) priors switchingrates where
                                                 mats = map (\x-> snd $ x [bl,sigma]) models
                                                 switchingrates = map (switchingSum (dataSize dt)) $ zip pis mats 

switchingSum nc (pi,mat) = getSwitchingRate mat pi nc

                                                                                
jcS :: Matrix Double
jcS = gtrS (replicate 5 1.0)

jcPi :: Vector Double
jcPi = fromList $ replicate 4 0.25

jcF = customF jcS jcPi

zeroParam v = (\x -> v,0,[],[])

genS a b c d e f = (4><4) [0,a,b,c,a,0,d,e,b,d,0,f,c,e,f,0]

gtrS :: [Double] -> Matrix Double
gtrS [a,b,c,d,e] = genS a b c d e 1
gtrPi :: [Double] -> Vector Double
gtrPi = piByLog

hkyS [a,b] = genS (a+b) b b b b (a+b)
hkyPi = piByLog


piByLog xs = fromList $ normaliseL (1.0:(map exp xs))
logByPi (first:xs) = let mlnFirst = -(log first) in
                        map (+mlnFirst) $ map log xs


type SPiFunctionTuple = (([Double] -> Matrix Double,Int,[Maybe Double],[Maybe Double]),([Double] -> Vector Double,Int,[Maybe Double],[Maybe Double]))
gtrF :: SPiFunctionTuple
gtrF = ((gtrS,5,replicate 5 $ Just 0.00001,replicate 5 Nothing),(gtrPi,3,replicate 3 $ Just (-20.0),replicate 3 $ Just 20.0))
hkyF :: SPiFunctionTuple
hkyF = ((hkyS,2,replicate 2 $ Just 0.00001, replicate 2 Nothing),(hkyPi,3,replicate 3 $ Just (-20.0),replicate 3 $ Just 20.0))

customF :: Matrix Double -> Vector Double -> SPiFunctionTuple
customF s pi = (zeroParam s,zeroParam pi)
wagF :: SPiFunctionTuple
wagF = customF wagS wagPi
jttF :: SPiFunctionTuple
jttF = customF jttS jttPi

getSensibleParams (f,n,(lower:lowers),(upper:uppers)) = (start:(getSensibleParams (f,n-1,lowers,uppers))) where
         start = case (lower,upper) of
                    (Nothing,Nothing) -> 0.0
                    (Just x,Just y) -> (x+y)/2
                    (Just x,Nothing) -> x+1.0
                    (Nothing,Just x) -> x-1.0
getSensibleParams (_,0,_,_)=[]
