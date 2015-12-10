import Text.ParserCombinators.Parsec
import System.IO
import System.Environment
import Data.List
import qualified Data.ByteString.Lazy as BL
import Data.Binary.Put
        
valouevFile :: GenParser Char st [(String, String, String, [Float])]
valouevFile =
            do  result <- many record
                eof
                return result

record :: GenParser Char st (String, String, String, [Float])
record = do map_name <- fieldContent ; eol
            enz_name <- parseEnzyme ; char '\t'; enz_acr_name <- fieldContent; char '\t'; frags <- sepBy fieldContent (char '\t'); eol
            eol
            return (map_name, enz_name, enz_acr_name, fmap read frags)

parseEnzyme :: GenParser Char st String
parseEnzyme = do (char '\t' >> fieldContent)
                 <|>
                 fieldContent
                 


fieldContent :: GenParser Char st String
fieldContent = many (noneOf "\t\n")


-- The end of line character is \n
eol :: GenParser Char st Char
eol = char '\n'
            
parseOM input = parse valouevFile "(unknown)" input



extract_frags :: (String, String, String, [Float]) -> [Float]
extract_frags (_, _, _, frags) = frags

frag_delim :: [Float]
frag_delim = [100000.0]

serialiseSomething :: [[Float]] -> Put
serialiseSomething skip_nodes = do
  let num_nodes = fromIntegral $ sum $ fmap length skip_nodes
  putWord32host num_nodes

             
main = do
  args <- getArgs
  let fname = (head args)
  let ofname = (last args)
  hdl <- openFile fname ReadMode 
  ohdl <- openFile ofname WriteMode
  contents <- hGetContents hdl
  case parseOM contents of
   Left x -> print $ show $ x
   Right x -> print $ show $ intercalate frag_delim $ fmap extract_frags x
  case parseOM contents of
   Left x -> print $ show $ x
   Right x -> print $ show $ intercalate frag_delim $ fmap extract_frags x
  BL.hPut ohdl $ runPut $ serialiseSomething [[1,2], [2,3]]
  hClose ohdl  

-- thanks to http://stackoverflow.com/questions/27726739/implementing-an-efficient-sliding-window-algorithm-in-haskell
windows :: Int -> [a] -> [[a]]
windows m = transpose . take m . tails

nth_skipnodes n nodes  = fmap sum $ windows n nodes

