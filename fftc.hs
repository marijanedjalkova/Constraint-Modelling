import Data.Complex
import Samples
import System.Environment
import Control.Parallel

-- twiddle factors
tw :: Int -> Int -> Complex Float
tw n k = cis (-2 * pi * fromIntegral k / fromIntegral n)

-- Discrete Fourier Transform -- O(n^2)
dft :: [Complex Float] -> [Complex Float]
dft xs = firstList `par` secondList `par` (firstList ++ secondList)
    where n = length xs
          firstList = [(calcDFT2 xs k) | k<-[0..div n 2] ]
          secondList = [(calcDFT2 xs k) | k<-[(div n 2) + 1..n-1] ]



calcDFT2::[Complex Float]-> Int -> Complex Float
calcDFT2 xs k = (subList k) `par` sum (subList k)
    where subList k = calcSubList xs (length xs) 0 k

{- This works but is not necessarily task control. Can argue that it is a case of divide and conquer though
dft :: [Complex Float] -> [Complex Float]
dft xs = calcDFT xs 0 

calcDFT::[Complex Float]-> Int -> [Complex Float]
calcDFT xs k 
    | k > n-1 = []
    | otherwise = newsum `par` rest `par` (newsum : rest)
    where n = length xs
          newsum = (subList k) `par` sum (subList k)
          subList k = calcSubList xs n 0 k
          rest = calcDFT xs (k+1)
-}
calcSubList :: [Complex Float] -> Int-> Int-> Int -> [Complex Float]
calcSubList [] _ _ _ = []
calcSubList [x] n j k = [x * (tw n (j*k))]
calcSubList (x:xs) n j k = h : rest
    where h = x * (tw n (j*k))
          rest = calcSubList xs n (j+1) k
-- Fast Fourier Transform

-- In case you are wondering, this is the Decimation in Frequency (DIF) 
-- radix 2 Cooley-Tukey FFT

fft :: [Complex Float] -> [Complex Float]
fft [a] = [a]
fft as = interleave ls rs
  where
    (cs,ds) = bflyS as
    ls = fft cs
    rs = fft ds

interleave [] bs = bs
interleave (a:as) bs = a : interleave bs as

bflyS :: [Complex Float] -> ([Complex Float], [Complex Float])
bflyS as = (los,rts)
  where
    (ls,rs) = halve as
    los = zipWith (+) ls rs
    ros = zipWith (-) ls rs
    rts = zipWith (*) ros [tw (length as) i | i <- [0..(length ros) - 1]]


-- split the input into two halves
halve as = splitAt n' as
  where
    n' = div (length as + 1) 2

-- the main function
-- uses Samples.samples to generate some sample data
--   samples :: Int -> Int -> [Complex Float]

defsize = 1000 -- change this to get larger samples
defseed = 1

main = do args <- getArgs
          let arglen = length args
          let n = argval args 0 defsize
          let seed = argval args 1 defseed
          let fun = if arglen > 2 && args !! 2 == "dft" then dft else fft
          print (sum (fun (samples seed n)))

argval args n def = if length args > n then
                       read (args !! n)
                     else def