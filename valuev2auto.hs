import Text.ParserCombinators.Parsec
import System.IO
import System.Environment
import Data.List

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

mylen :: [a] -> Int
mylen [] = 0
mylen (x:xs) = 1 + mylen xs

joinem :: [(String, String, String, [Float])] -> [Float]
joinem [] = []
joinem (r:rs) = let (name, enz_name, enz_acr, frags) = r
              in frags ++ [1000000] ++ (joinem rs) -- let (name, enz_name, enz_acr, frags) = r
--
extract_frags :: (String, String, String, [Float]) -> [Float]
extract_frags (_, _, _, frags) = frags

frag_delim :: [Float]
frag_delim = [100000.0]

main = do
  args <- getArgs
  let fname = (head args)
  hdl <- openFile fname ReadMode 
  contents <- hGetContents hdl
  case parseOM contents of
   Left x -> print $ show $ x
   Right x -> print $ show $ intercalate frag_delim $ fmap extract_frags x
  
