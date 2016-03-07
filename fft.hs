import Data.Complex
import Samples
import System.Environment
import Control.Parallel

-- twiddle factors
tw :: Int -> Int -> Complex Float
tw n k = cis (-2 * pi * fromIntegral k / fromIntegral n)

-- Discrete Fourier Transform -- O(n^2)
dft :: [Complex Float] -> [Complex Float]
dft xs = [ sum [ xs!!j * tw n (j*k) | j <- [0..n']] | k <- [0..n']]
  where
    n = length xs
    n' = n-1


--parZipWith ::  (a -> b -> c) -> [a] -> [b] -> [c]
--parZipWith strat z as bs = 
 -- zipWith z as bs `using` parList strat

--parZipWith :: ( Traversable t , Foldable f , NFData c ) => ( a −> b −> c ) −> t a −> f b −> Par ( t c )
--parZipWith comb t f = let step ( x : xs ) y = ( xs , spawnP $ comb y x ) in 
   -- (T.sequence . snd $ T.mapAccumL step (F.toList f ) t ) >>= T.mapM get



-- Fast Fourier Transform

-- In case you are wondering, this is the Decimation in Frequency (DIF) 
-- radix 2 Cooley-Tukey FFT

fft :: [Complex Float] -> [Complex Float]
fft [a] = [a]
fft as = my_interleave ls rs
  where
    (cs,ds) = bflyS as
    ls = fft cs
    rs = fft ds

my_interleave [] bs = bs
my_interleave bs [] = bs
my_interleave (a:as) (b:bs) = let h = [a] ++ [b] in
    h `par` (my_interleave as bs)


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

-- 