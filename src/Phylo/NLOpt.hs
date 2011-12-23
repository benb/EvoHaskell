{-# LANGUAGE ForeignFunctionInterface #-}
module Phylo.NLOpt (bobyqa,cobyla) where

import Data.Maybe (fromMaybe)
import Foreign
import Foreign.C
import Foreign.C.Types
import Foreign.Ptr (Ptr, FunPtr, freeHaskellFunPtr)
import Control.Monad

foreign import ccall safe "nlopt_c opt_bobyqa" 
        bobyqa_ :: CDouble -> Ptr CDouble -> CUInt -> FunPtr ((CUInt -> Ptr CDouble -> Ptr CDouble -> Ptr () -> IO CDouble)) -> Ptr CDouble -> Ptr CDouble -> CInt
foreign import ccall safe "nlopt_c opt_cobyla" 
        cobyla_ :: CDouble -> Ptr CDouble -> CUInt -> FunPtr ((CUInt -> Ptr CDouble -> Ptr CDouble -> Ptr () -> IO CDouble)) -> Ptr CDouble -> Ptr CDouble -> CInt

foreign import ccall "wrapper"
        wrap :: (CUInt -> Ptr CDouble -> Ptr CDouble -> Ptr () -> IO CDouble) -> IO (FunPtr ((CUInt -> Ptr CDouble -> Ptr CDouble -> Ptr () -> IO CDouble)))

nlopt :: (CDouble -> Ptr CDouble -> CUInt -> FunPtr ((CUInt -> Ptr CDouble -> Ptr CDouble -> Ptr () -> IO CDouble)) -> Ptr CDouble -> Ptr CDouble -> CInt) -> Double -> [Double] -> ([Double] -> Double) -> [Maybe Double] -> [Maybe Double] -> IO ([Double],Int)

bobyqa = nlopt bobyqa_
cobyla = nlopt cobyla_
nlopt met xtol params f lower upper = do lower' <- newArray $ map (realToFrac . fromMaybe (-1E100)) lower
                                         upper' <- newArray $ map (realToFrac . fromMaybe 1E100) upper 
                                         let np = length params
                                         let f' a b c d = fmap realToFrac $ fmap f (fmap (map realToFrac) $ peekArray np b)
                                         f'' <- wrap f'
                                         startP <- newArray $ map (realToFrac) params
                                         let retCode = fromIntegral $ met (realToFrac xtol) startP (fromIntegral np) f'' lower' upper' 
                                         ans <- seq retCode $ fmap (map realToFrac) $ peekArray np startP
                                         return (ans,retCode)

