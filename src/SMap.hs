import System.Environment (getArgs)
import System.Console.GetOpt
import Phylo.Alignment
import Data.List
import Control.Monad
import Phylo.Tree
import Phylo.Data
import Phylo.Likelihood
import Phylo.Opt
import Numeric.GSL.Minimization
import Data.Packed.Vector
import Debug.Trace
import System.Random
import Data.Packed.Matrix
import System.IO
import System.Posix.Types
import System.Posix.IO
import Foreign.C.String (CString,newCString)
import Phylo.Matrix
import Stochmap
import Data.Char (isSpace)
import Phylo.NeXML
import Statistics.Quantile
import Phylo.Graphics.Plotting
import Graphics.Rendering.Chart.Renderable (renderableToPDFFile)
import qualified Data.Vector.Unboxed as UVec

data Flag = NumCats String | AlignmentFile String | TreeFile String | Inter | Intra 
        deriving (Show,Eq,Ord)
options = [ Option ['a'] ["alignment"] (ReqArg AlignmentFile "FILE") "Alignment",
            Option ['t'] ["tree"] (ReqArg TreeFile "FILE") "Tree",
            Option [] ["num-cats"] (ReqArg NumCats "NUMGAMMACATS") "Number of Gamma Rate Categories" ,
            Option [] ["inter"] (NoArg Inter) "Inter",
            Option [] ["intra"] (NoArg Intra) "Intra"]

main = do args <- getArgs
          (aln,tree,stdGen,remainOpts,nonOpts) <- case getOpt Permute options args of 
                         (((AlignmentFile aln):(TreeFile tre):xs),nonOpts,[]) -> do aln <- parseAlignmentFile parseUniversal aln
                                                                                    tree <- (liftM readBiNewickTree) (readFile tre)
                                                                                    stdGen <- getStdGen
                                                                                    return (aln,tree,stdGen,xs,nonOpts)
                         (_,_,msgs) -> error $ concat msgs
          let (cats,ii,remainOpts') = case remainOpts of 
                (NumCats n):i:rem -> (read n,i,rem)
                (i:rem) -> (4,i,rem)
          let lMat = case ii of
                Inter -> interLMat (cats+1) 20
                Intra -> intraLMat (cats+1) 20
          output <- case (aln,tree,nonOpts) of 
                (Just a,Right t,params)->   do let alpha:sigma:priorZero:[] = map read $ take 3 params
                                               let pAln = pAlignment a
                                               let (nSite,nCols,multiplicites) = getAlnData pAln where
                                                   getAlnData (PatternAlignment names seqs columns patterns counts) = (length counts, length columns,counts)
                                               let piF = fromList $ safeScaledAAFrequencies a
                                               let model = thmmModel (cats+1) wagS piF [priorZero,alpha,sigma]
                                               let t2 = addModelFx (structDataN (cats+1) AminoAcid (pAln) t) model [1.0]
                                               let nState = (cats+1)*20
                                               let nProc = 1 {-- fixme --}
                                               let pi_i = map toList (snd model) 
                                               let mixProbs = getPriors t2
                                               let qMat = thmmModelQ (cats+1) wagS piF [priorZero,alpha,sigma]
                                               let qset = [toLists qMat]
                                               let pBEStr = getPartialBranchEnds t2
                                               let partials = map (map (map (toLists . trans))) $ toPBEList pBEStr
                                               let branchLengths = toPBELengths pBEStr
                                               let nBranch =  length $ toPBELengths pBEStr
                                               let sitelikes = [likelihoods t2] {-- FIXME, outer dimension is nProc, assume 1 --}
                                               let sitemap = mapBack pAln
                                               putStrLn $ show $ zip [0..] $ getCols pAln 
                                               let xml = xmlTree t2
                                               putStrLn "Here comes the xml:"
                                               putStr xml
                                               let stochmapF a b c d e f g = calculateAndWrite nSite nState nBranch nProc nCols lMat multiplicites a b c d e f g Nothing
                                               let stochmapT tree = stochmapF sitemap' partials' qset' sitelikes' pi_i' branchLengths' mixProbs' where
                                                                    sitemap' = mapBack $ pAlignment pAln'
                                                                    pAln' = getAln tree
                                                                    pBEStr' = getPartialBranchEnds tree
                                                                    partials' = map (map (map (toLists . trans))) $ toPBEList pBEStr'
                                                                    qset' = map toLists $ getQMat tree
                                                                    sitelikes' = [likelihoods tree] --FIXME
                                                                    pi_i' = map toList $ getPi tree
                                                                    branchLengths' = toPBELengths pBEStr
                                                                    mixProbs' = getPriors tree
                                               ans <- stochmapT t2 -- sitemap partials qset sitelikes pi_i branchLengths mixProbs Nothing
                                               putStrLn "OK0"
                                               let numSim = 10
                                               let numQuantile = 100
                                               let stdGens = take numSim $ genList stdGen
                                               print stdGens
                                               let alnLength = length $ Phylo.Likelihood.columns pAln
                                               putStrLn "OK0.1"
                                               print $ map next stdGens
                                               let simulations = map optBLDFull0 $ map (\x->makeSimulatedTree AminoAcid (cats+1) x alnLength t2) stdGens
                                              -- let simulations =  map (\x->t2) stdGens
                                               let simPriors = map getPriors simulations
                                               putStrLn $ "OK0.2 "  ++ (show $ length simulations)
                                               print $ simulations
                                               print $ map logLikelihood simulations
                                               simStoch <- mapM stochmapT simulations
                                               let simMappings = map treeToMap simulations
                                               putStrLn $ "OK0.3 "  ++ (show $ length simMappings)
                                               putStrLn $ "OK0.4 "  ++ (show $ length simStoch)
                                               putStrLn "OK1"
                                               print $ drop 9 simStoch
                                               putStrLn "OK2"
                                               putStrLn $ "sMappings " ++ (show $ length simMappings)
                                               putStrLn $ "sMappings " ++ (show $ length simPriors)
                                               putStrLn $ "sMappings " ++ (show $ length simStoch)
                                               let quantileDist = qList (uniformQRaw simMappings simPriors simStoch) numQuantile
                                               print quantileDist 
                                               let revQuantile = revQuantileF ans numQuantile quantileDist t2
                                               let simQuantiles = map (\(x,y) ->  revQuantileF x numQuantile quantileDist y ) $ zip simStoch simulations
                                               let ninetyFive = map (\x-> (x!!a,x!!b)) $ map sort $ transpose simQuantiles where
                                                            b :: Int
                                                            a :: Int
                                                            a = floor $ (fromIntegral numSim) * 0.025
                                                            b = (ceiling ((fromIntegral numSim) * 0.975)) - 1
                                               putStrLn $ "uniform? " ++ (show revQuantile)
                                               renderableToPDFFile (makePlot revQuantile ninetyFive PDF) 640 480 "out.pdf"
                                               --stochmapOut ans sitemap [1.0] putStr
                                               print "OK"
                                               return Nothing
          case output of
               Just str -> putStrLn str
               Nothing -> return ()

revQuantileF x numQuantile quantileDist t2 = (flip qList) numQuantile $ myQuantile $  map (/(fromIntegral numQuantile)) $ map fromIntegral $ map (reverseQuantile quantileDist) $ linearFromRaw (mapBack $ pAlignment $ getAln t2) (getPriors t2) x
treeToMap tree = mapBack $ pAlignment $ getAln tree

allInGeneric method sets (l,r) = case (findIndex (==True) $ map (method l r) sets) of
                                        Just a -> a
                                        Nothing -> other
                                 where other = length sets
allIn = allInGeneric allInDynamic'
allInNoRoot = allInGeneric allInNoRoot'


allIn' :: [String] -> [String] -> [String] -> Bool
allIn' l r set = (l \\ set == [] || r \\ set == []) --we want l or r to be fully contained within 'set'

allInNoRoot' :: [String] -> [String] -> [String] -> Bool
allInNoRoot' l r set = ((l \\ set == []) && (set \\ l /= [])) || ((r \\ set == []) &&  (set \\ r /= [])) --we want l or r to be fully contained with in 'set', but not the same as 'set'
                        
--if first symbol is '@' then use noRoot version
allInDynamic' :: [String] -> [String] -> [String] -> Bool
allInDynamic' l r ("@":set) = allInNoRoot' l r set
allInDynamic' l r set = allIn' l r set


splitBy delim s = ans where
  (token,rest) = span (/=delim) s
  ans = case rest of
    [] -> [token]
    t  -> token : splitBy delim (tail rest)

splitsStr = unlines . splitsStr' 0
splitsStr' i [] = []
splitsStr' i ((l,r):xs) = ((show i) ++ " " ++ (intercalate " " l) ++ " | " ++ (intercalate " " r)) : (splitsStr' (i+1) xs)

getCols (PatternAlignment _ _ c _ _) = c

normalise list = map ( / total) list where     
                 total = foldr (+) 0.0 list

safeScaledAAFrequencies = normalise . map (\x-> if x < 1e-15 then 1e-15 else x) . scaledAAFrequencies
optParamsAndBL model tree params priors lower upper cutoff = optWithBS' [] cutoff (0,0) Nothing lower upper priors model (dummyTree tree) params                                                            
optParamsAndBLIO model tree params priors lower upper cutoff = optWithBSIO' [] cutoff (0,0) Nothing (map (\x->0.01) lower) (map (\x->1E-4) lower) lower upper priors model (dummyTree tree) params                                                            

trim = f . f where 
   f = reverse . dropWhile isSpace
clean = clean' . trim where
   clean' (' ':' ':xs) = clean' (' ':xs)
   clean' (x:xs) = x:(clean' xs)
   clean' [] = []
joinWith i (x:x':xs) = (show x) ++ i ++ (joinWith i (x':xs))
joinWith i [x] = (show x) 
joinWith i [] = []

stochmapOrder :: ([[[Double]]],[[Double]]) -> [Int] -> [Double] -> ([[Double]],[Double])
stochmapOrder (condE,priorE) mapping priors = order where
                                                condE' = map (fixProc2 priors) condE
                                                priorE' = map (fixProc priors) priorE
                                                fixProc pr x = (foldr (+) (0.0) (map (\(x,y) -> x*y) $ zip pr x))
                                                fixProc2 pr xs = (map $ fixProc pr) (transpose xs)
                                                order = (transpose $ map (\i-> (map (!!i) condE')) mapping, priorE')

stochmapOut :: ([[[Double]]],[[Double]]) -> [Int] -> [Double] -> (String -> IO()) -> IO ()
stochmapOut (condE',priorE') mapping priors f = do 
                                                f header
                                                (mapM . mapM) f remainder
                                                f "\n"
                                                f headerSite 
                                                mapM f remainderSite
                                                f "\n"
                                                f headerBranch
                                                mapM f remainderBranch
                                                return () where
                                                 (condE''',priorE''') = stochmapOrder (condE',priorE') mapping priors
                                                 priorE'' = replicate (length mapping) priorE'''
                                                 condE'' = transpose condE'''
                                                 condE=transpose condE''
                                                 priorE=transpose priorE''
                                                 header = "Branch\tSite\tConditional_expectation\tPrior_expectation\n"
                                                 remainder = fmtBranch [0..] $ zip condE priorE 
                                                 fmtBranch (b:bs) ((cond,prior):xs) = (fmtBranch' b [0..] cond prior)  : (fmtBranch bs xs)
                                                 fmtBranch bs [] = []
                                                 fmtBranch' b sites [] priors = []
                                                 fmtBranch' b (site:sites) (cond:cs) (prior:ps) = ((show b) ++ "\t" ++ (show site) ++ "\t" ++ (show cond) ++ "\t" ++ (show prior) ++ "\n") : (fmtBranch' b sites cs ps)
                                                 headerBranch = "Branch\tTotal_conditional_expectation\tTotal_prior_expectation\n"
                                                 remainderBranch = fmtBranch2 [0..] $ zip condE priorE
                                                 total xs = foldr (+) 0.0 xs
                                                 fmtBranch2 bs [] = []
                                                 fmtBranch2 (b:bs) ((cond,prior):xs) = ((show b) ++"\t" ++ (show $ total cond) ++ "\t" ++ (show $ total prior) ++ "\n") : (fmtBranch2 bs xs)
                                                 headerSite = "Site\tTotal_conditional_expectation\tTotal_prior_expectation\n"
                                                 remainderSite = fmtBranch2 [0..] $ zip condE'' priorE''

qList qFunc numQ = map (\x-> qFunc x numQ) [0..numQ]

myQuantile myData q k = continuousBy medianUnbiased q k (UVec.fromList myData)


reverseQuantile myQList point = fst $ head $ filter (\(x,y) -> point <= y) $ zip [0..] myQList 

uniformQ :: [([[Double]],[Double])] -> (Int -> Int -> Double)
uniformQ simData = partialFunc where 
                      partialFunc q k | trace ("qData " ++ (show q ) ++ " " ++ (show k) ++ " : " ++  (show $ concatMap allDisc simData)) True = continuousBy medianUnbiased q k (UVec.fromList $ concatMap allDisc simData)

linearFromRaw mapping priors ans = linearAns $ stochmapOrder ans mapping priors
linearAns reformattedAns = allDisc reformattedAns
reformatAns mapping priors stochResult = (map $ uncurry3 stochmapOrder) $ zip3 stochResult mapping priors
uniformQRaw mapping priors stochResult = uniformQ $ reformatAns mapping priors stochResult

uncurry3 f (a,b,c) = f a b c 

allDisc :: ([[Double]],[Double]) -> [Double]
allDisc (condE:xs,priorE:ys) | traceShow ("len " ++ (show $ 1 + (length xs)) ++ " " ++ (show $ 1 + (length ys))) True = (allDisc' condE priorE) ++ (allDisc (xs,ys))
allDisc ([],[]) = []
allDisc (c,p) | traceShow (c,p) True = undefined
allDisc' condE priorE | trace ("prior cond " ++ (show priorE) ++ "  " ++ (show condE)) True = map (\x -> priorE - x) condE

