import Data.Complex
import Samples
import System.Environment
import Control.Parallel

-- twiddle factors
tw :: Int -> Int -> Complex Float
tw n k = cis (-2 * pi * fromIntegral k / fromIntegral n)

-- Discrete Fourier Transform -- O(n^2)
dft :: [Complex Float] -> [Complex Float]
dft xs = n`par` n' `par` [ sum [ xs!!j * tw n (j*k) | j <- [0..n']] | k <- [0..n']]
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
-- 1. bflyS as - divide in half, zip sum and diff into two lists. The first one ready to return, the second one evald with tws and zipped and returned
-- 2. ls = fft cs
-- 3. rs = fft ds
-- (2 and 3 can be parallelised) - including 1!!!
-- 4. interleave - yes, should stay parallelised, probably.
-- basically can do everything together, getting and interleaving values as they appear
fft :: [Complex Float] -> [Complex Float]
fft [a] = [a]
fft as = (cs,ds) `par` ls `par` rs `par` par_interleave ls rs
  where
    (cs,ds) = bflyS as
    ls = fft cs
    rs = fft ds

par_interleave [] bs = bs
par_interleave bs [] = bs
par_interleave (a:as) (b:bs) = h `par` (rest `pseq` (h ++ rest))
    where h = [a,b]
          rest = par_interleave as bs

-- 1. halve as - no need to parallelise
-- 2. los
-- 3. ros (los and ros can be task parallelised)
-- 4. rts - list comprehension can be parallelised
bflyS :: [Complex Float] -> ([Complex Float], [Complex Float])
bflyS as = los `par` ros `par` (los,rts)
  where
    (ls,rs) = halve as
    los = parZipWith (+) ls rs
    ros = parZipWith (-) ls rs
    list = [0..(length ros) - 1]
    tweaked = tweakList (length as) list
    rts = parZipWith (*) ros tweaked 

tweakList :: Int -> [Int] -> [Complex Float]
tweakList length [] = []
tweakList length [x] = [tw length x]
tweakList length (x:xs) = h `par` (rest `pseq` (h : rest))
    where h = tw length x
          rest = tweakList length xs

parZipWith :: (a->a->a) -> [a] -> [a] -> [a]
parZipWith f [] _ = []
parZipWith f [x] list = zipWith f [x] list
parZipWith f (x:xs) (y:ys) = h `par` (rest `pseq` (h : rest))
    where h = f x y
          rest = parZipWith f xs ys

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