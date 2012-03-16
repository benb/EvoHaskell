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

data Flag = NumCats String | AlignmentFile String | TreeFile String | Alpha String | OptAlpha | OptThmm | OptThmmP | OptThmm2 String | OptSim1 | OptSim2 String | OptSim0 String | Seed String | ThmmStochmap String | NoOpt 
        deriving (Show,Eq,Ord)
options = [ Option ['a'] ["alignment"] (ReqArg AlignmentFile "FILE") "Alignment",
            Option ['t'] ["tree"] (ReqArg TreeFile "FILE") "Tree",
            Option ['s'] ["seed"] (ReqArg Seed "INTEGER") "Seed",
            Option [] ["stochmap"] (ReqArg ThmmStochmap "PARAMS") "Calculate stochastic mapping for Thmm",
            Option [] ["num-cats"] (ReqArg NumCats "NUMGAMMACATS") "Number of Gamma Rate Categories",
            Option ['g'] ["gamma"] (ReqArg Alpha "DOUBLE") "Use Gamma Model" ,
            Option ['o'] ["opt-gamma"] (NoArg OptAlpha) "Optimise Gamma Model",
            Option ['m'] ["opt-thmm"] (NoArg OptThmm) "Optimise THMM Model" ,
            Option ['n'] ["opt-thmmplus"] (NoArg OptThmmP) "Optimise THMM Model",
            Option ['p'] ["opt-thmm2"] (ReqArg OptThmm2 "splitsStr") "2 state THMM",
            Option [] ["noopt"] (NoArg NoOpt) "Don't optimize",
            Option [] ["sim0"] (ReqArg OptSim0 "PARAMS") "Simulate WAG",
            Option [] ["sim1"] (NoArg OptSim1) "Simulation 1",
            Option [] ["sim2"] (ReqArg OptSim2 "PARAMS") "Simulation 2"]

main = do args <- getArgs
          (aln,tree,stdGen,remainOpts) <- case getOpt Permute options args of 
                         (((AlignmentFile aln):(TreeFile tre):(Seed seed):xs),[],[]) -> do aln <- parseAlignmentFile parseUniversal aln
                                                                                           tree <- (liftM readBiNewickTree) (readFile tre)
                                                                                           let stdGen = mkStdGen $ read seed
                                                                                           return (aln,tree,stdGen,xs)
                         (((AlignmentFile aln):(TreeFile tre):xs),[],[]) -> do aln <- parseAlignmentFile parseUniversal aln
                                                                               tree <- (liftM readBiNewickTree) (readFile tre)
                                                                               stdGen <- getStdGen
                                                                               return (aln,tree,stdGen,xs)
                         (((TreeFile tre):(Seed seed):xs),[],[]) -> do tree <- (liftM readBiNewickTree) (readFile tre)
                                                                       let stdGen = mkStdGen $ read seed
                                                                       return (Nothing,tree,stdGen,xs)
                         (_,_,msgs) -> error $ concat msgs
          let (cats,remainOpts') = case remainOpts of 
                (NumCats n):rem -> (read n,rem)
                rem -> (4,rem)
          output <- case (aln,tree,(sort remainOpts')) of 
                (Just a,Right t,[])-> return $ Just $ "WAG lkl:" ++ (show $ quickLkl a t wagPi wagS)
                (Just a,Right t,(Alpha val):[])-> return $ Just $ "WAG Gamma lkl:" ++ (show $ quickGamma cats (read val) a t wagPi wagS)
                (Just a,Right t,(OptAlpha):[])-> return $ Just $ "Opt Alpha: " ++ (show alpha) ++ " " ++ (show lkl) where
                                                        t2 = structDataN 1 AminoAcid (pAlignment a) t
                                                        piF = fromList $ safeScaledAAFrequencies a
                                                        model = gammaModel cats wagS piF
                                                        (optTree,[logalpha]) = optParamsAndBL model t2 [0.5] [0.25,0.25,0.25,0.25] [Just 0.001] [Nothing] 0.01
                                                        alpha = exp logalpha
                                                        lkl = logLikelihood optTree
                (Just a,Right t,(OptThmm):[])-> do 
                        let t2 = structDataN (cats+1) AminoAcid (pAlignment a) t                                                                                                                                                                       
                        let piF = fromList $ normalise $ map (\x-> if x < 1e-15 then 1e-15 else x) $ safeScaledAAFrequencies a
                        putStrLn "OK"
                        let model = thmmModel (cats+1) wagS piF
                        (optTree',[priorZero,alpha,sigma]) <- optParamsAndBLIO model t2 [0.21546728749044514,0.8209232358343468,10.33751049838077] [1.0] [Just 0.01,Just 0.001, Just 0.00] [Just 0.99,Just 100.0,Just 10000.0] 0.01
                        let optTree = annotateTreeWithNumberSwitchesSigma sigma optTree'
                        let lkl = logLikelihood optTree
                        return $ Just $ "Opt Thmm: " ++ (show alpha) ++ " " ++ (show sigma) ++ " " ++ (show priorZero) ++ " " ++ (show lkl) ++ "\n" ++ (show optTree) 
                (Just a,Right t,(ThmmStochmap params ):[])-> do let alpha:sigma:priorZero:[] = map read $ take 3 $ words params
                                                                let interintra = head $ drop 3 $ words params
                                                                let pAln = pAlignment a
                                                                let (nSite,nCols,multiplicites) = getAlnData pAln where
                                                                    getAlnData (PatternAlignment names seqs columns patterns counts) = (length counts, length columns,counts)
                                                                let piF = fromList $ safeScaledAAFrequencies a
                                                                let model = thmmModel (cats+1) wagS piF [priorZero,alpha,sigma]
                                                                let t2 = addModelFx (structDataN (cats+1) AminoAcid (pAln) t) model [1.0]
                                                                let nState = (cats+1)*20
                                                                let lMat = case interintra of 
                                                                                "inter" -> interLMat (cats+1) 20
                                                                                "intra" -> intraLMat (cats+1) 20
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
                                                                --handle <- openFile "out.txt" WriteMode
                                                                putStrLn $ show $ zip [0..] $ getCols pAln 
                                                                putStr $ splitsStr $ toPBESplits pBEStr
                                                                --putStrLn $ "OK " ++ (show $ accEigQ (qMat) ((snd model)!!0))
                                                                calculateAndWrite nSite nState nBranch nProc nCols lMat multiplicites sitemap partials qset sitelikes pi_i branchLengths mixProbs stdout
                                                                return Nothing
                                                                {-c_test chandle $ fromIntegral 10-}
                (Just a,Right t,(OptThmmP):[])-> return $ Just $ "Opt Thmm: " ++ (show alpha) ++ " " ++ (show optTree) ++ " " ++ (show priorZero) ++ " " ++ (show lkl) where
                                                        piF = fromList $ safeScaledAAFrequencies a
                                                        t2 = addModelFx (structDataN (cats+1) AminoAcid (pAlignment a) t) (gammaModel cats wagS piF [0.5]) $ flatPriors cats                                                                                                                                                                
                                                        startBL = getBL t2
                                                        newBL = map (\x -> x:0.0:[]) startBL
                                                        t3 = fst $ setBLX' 0 newBL t2
                                                        model = thmmPerBranchModel (cats+1) wagS piF
                                                        (optTree,[priorZero,alpha]) = optParamsAndBL model t3 [0.1,0.5] [1.0] [Just 0.01,Just 0.001] [Just 0.99,Nothing] 0.01
                                                        bls = getBL t3
                                                        lkl = logLikelihood optTree
                (Just a,Right t,(OptThmm2 spl'):xs)-> do (optTree',optParams) <- case xs of 
                                                                list | NoOpt `elem` list -> return (t4,p) where
                                                                                              p = (initSigma ++ [initP0,initAlpha])
                                                                                              t4 | trace ("NUMMODELS " ++ (show numModels)) True = getFuncT1A [1.0] mapped (1,numModels) model t3 p
                                                                _ -> optBSParamsBLIO (1,numModels) mapped (map (\x->0.01) lower) lower upper [1.0] model t3 (initSigma ++ [initP0,initAlpha])
                                                         let optTree = annotateTreeWithNumberSwitches optTree'
                                                         print optTree
                                                         let (sigma,[priorZero,alpha]) = splitAt numModels optParams
                                                         let lkl = logLikelihood optTree
                                                         putStrLn ("SIGMA " ++ (show sigma))
                                                         --
                                                         let (nSite,nCols,multiplicites) = getAlnData (pAlignment a) where
                                                                    getAlnData (PatternAlignment names seqs columns patterns counts) = (length counts, length columns,counts)
                                                         let nState = (cats+1)*20
                                                         let pi_i = map toList $ getpi optTree' where
                                                                        getpi (DTree _ _ _ _ _ _ pi) = pi
                                                         let sitelikes = [likelihoods optTree']
                                                         let nProc = 1 {-- fixme --}
                                                         let pBEStr = getPartialBranchEnds optTree'
                                                         let qset = map toLists $ map head $ toPBEQ pBEStr
                                                         let mixProbs = getPriors optTree'
                                                         let nBranch =  length $ toPBELengths pBEStr
                                                         let branchLengths = toPBELengths pBEStr
                                                         let branchLengths2 = everyOther $ getBL optTree' where 
                                                                everyOther (x:y:z) = x:(everyOther z)
                                                                everyOther [] = []
                                                         print branchLengths 
                                                         print branchLengths2
                                                         --TODO TEST
                                                         let interintra = "inter"
                                                         let sitemap = mapBack $ pAlignment a
                                                         let partials = map (map (map (toLists . trans))) $ toPBEList pBEStr
                                                         let lMat = case interintra of 
                                                                   "inter" -> interLMat (cats+1) 20
                                                                   "intra" -> intraLMat (cats+1) 20
                                                         putStrLn $  "Opt Thmm: " ++ (show alpha) ++ " " ++ (joinWith " " $ sigma) ++ " " ++ " " ++ (show priorZero) ++ " " ++ (show lkl) ++ "\n" ++ (show optTree) 
                                                         putStrLn $ posteriorTipsCSV (cats+1) (pAlignment a) optTree'
                                                         --calculateAndWrite nSite nState nBranch nProc nCols lMat multiplicites sitemap partials qset sitelikes pi_i branchLengths mixProbs stdout
                                                         return Nothing where
                                                                spl = clean spl'
                                                                piF = fromList $ safeScaledAAFrequencies a
                                                                t2 = addModelFx (structDataN (cats+1) AminoAcid (pAlignment a) t) (gammaModel cats wagS piF [0.5]) $ flatPriors cats                                                                                                                                                             
                                                                startBL = getBL t2
                                                                newBL = map (\x -> x:0.0:[]) startBL
                                                                t3 = fst $ setBLX' 0 newBL t2
                                                                goodNodes = map (splitBy ',') $ drop 1 $ splitBy ' ' spl
                                                                mapped = makeMapping (allIn goodNodes) t3
                                                                numModels = length $ nub $ whichModels t3 mapped
                                                                initParams = map read $ splitBy ',' $ head $ splitBy ' ' spl
                                                                initAlpha = head initParams
                                                                initSigma = take numModels $ tail initParams
                                                                initP0 = last initParams
                                                                model | traceShow t3 True = thmmPerBranchModel (cats+1) wagS piF 
                                                                lower = (replicate numModels $ Just 0.0) ++ [Just 0.001,Just 0.001] 
                                                                upper = (replicate numModels Nothing) ++ [Just 0.99,Nothing]
                                                        {-optBSParamsBL numBSParam mapping = optWithBS' [] 1E-4 numBSParam (Just mapping)-}
                                                        {-optParamsBL = optWithBS' [] 1E-4 (0,0) Nothing-}
                                                        {-optWithBS' :: [(Double,IterType)] -> Double -> (Int,Int) -> (Maybe (DNode -> Int)) -> [Maybe Double] -> [Maybe Double] -> [Double] -> ModelF -> DNode -> [Double] -> (DNode,[Double])-}
                                                        {-optWithBS' iterations cutoff numBSParam mapping lower upper priors model tree startParams = (bestTree,bestParams) where-}
{-                (Just a,Right t,OptSim1:[]) -> putStrLn $ concat $ toFasta $ simulateSequences AminoAcid (cats+1) stdGen 349 t2 where-}
                                                        {-alpha = 0.8335109218715334-}
                                                        {-priorZero = 0.1907077998691572-}
                                                        {-sigma0 = 15.631076379213784-}
                                                        {-sigma1 = 1.0836675920110694e-8-}
                                                        {-piF = fromList $ scaledAAFrequencies a-}
                                                        {-t2 = addModelFx (setBLMapped 1 (dummyTree (structDataN (cats+1) AminoAcid (pAlignment a) t)) mapped ) (thmmPerBranchModel (cats+1) cpRevS cpRevPi [priorZero,alpha]) [1.0]-}
                                                        {-goodNodes = ["E_Nosloc","E_Enccun","E_Gluple","A_Aerper","A_Metbar"] -}
                                                        {-mapped = makeMapping (\(x,y) -> if (((x \\ goodNodes) == []) || ([] == (y \\ goodNodes))) then ([sigma0]) else ([sigma1])) t2-}
                (a,Right t,(OptSim0 params):[]) -> return $ Just $ concat $ toFasta $ simulateSequences AminoAcid (cats+1) stdGen (floor length) t2 where
                                                        alpha:sigma0:sigma1:priorZero:length:[] = map (read) $take 5 $ words params
                                                        goodNodes = splitBy ',' $ head $ drop 5 $ words params 
                                                        piF = case a of 
                                                                Just aln -> fromList $ safeScaledAAFrequencies aln
                                                                Nothing -> wagPi
                                                     --   piF | trace ((show sigma0) ++ "  " ++  (show sigma1)) True= fromList $ scaledAAFrequencies a
                                                        emptyAlignment = ListAlignment [] [] []
                                                        t2 = trace (show goodNodes) $ addModelFx (setBLMapped 1 (dummyTree (structDataN (cats+1) AminoAcid (pAlignment emptyAlignment) t)) mapped ) (thmmPerBranchModel (cats+1) wagS piF [priorZero,alpha]) [1.0]
                                                        mapped = makeMapping (\(x,y) -> if (((x \\ goodNodes) == []) || ([] == (y \\ goodNodes))) then ([sigma0]) else ([sigma1])) t2
                (Just a,Right t,OptSim2 params:[]) -> do (lgInnerS,lgInnerPi) <- parsePamlDatIO "lgInner"
                                                         (lgOuterS,lgOuterPi) <- parsePamlDatIO "lgOuter"
                                                         return $ Just $ compute (lgInnerS,lgInnerPi) (lgOuterS,lgOuterPi) where
                                                             compute (s0,pi0) (s1,pi1) = concat $ toFasta $ simulateSequences AminoAcid (cats+1) stdGen 349 t3 where
                                                                alpha:sigma0:sigma1:priorZero:[] = map (read) $ words params
                                                                piF = fromList $ safeScaledAAFrequencies a
                                                                (model0,pi_0) = thmmPerBranchModel (cats+1) s0 pi0 [priorZero,alpha]
                                                                (model1,pi_1) = thmmPerBranchModel (cats+1) s1 pi1 [priorZero,alpha]
                                                                t2 = addModelFx (setBLMapped 1 (dummyTree (structDataN (cats+1) AminoAcid (pAlignment a) t)) mapped ) (thmmPerBranchModel (cats+1) cpRevS cpRevPi [priorZero,alpha]) [1.0]
                                                                tfFunc x y = ((x \\ goodNodes) == []) || ([] == (y \\ goodNodes))
                                                                goodNodes = ["E_Nosloc","E_Enccun","E_Gluple","A_Aerper","A_Metbar"] 
                                                                mapped = makeMapping (\(x,y) -> if (tfFunc x y)  then ([sigma0]) else ([sigma1])) t2
                                                                mappedModels = makeMapping (\(x,y) -> if (tfFunc x y) then model1 else model0) t2
                                                                t3 = restructDataMapped t2 mappedModels [1.0] pi_0
                (_,_,_) -> error "Can't parse something"
          case output of
               Just str -> putStrLn str
               Nothing -> return ()


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
