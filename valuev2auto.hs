import Text.ParserCombinators.Parsec as PS
import System.IO
import System.Environment
import Data.List
import qualified Data.ByteString.Lazy as BL
import Data.Binary.Put
import Control.Applicative as A
    
valouevFile :: GenParser Char st [(String, String, String, [Float])]
valouevFile =
            do  result <- PS.many record
                eof
                return result

record :: GenParser Char st (String, String, String, [Float])
record = do map_name <- fieldContent ; eol
            enz_name <- parseEnzyme ; char '\t'; enz_acr_name <- fieldContent; char '\t'; frags <- sepBy fieldContent (char '\t'); eol
            eol
            return (map_name, enz_name, enz_acr_name, fmap read frags)

parseEnzyme :: GenParser Char st String
parseEnzyme = do (char '\t' >> fieldContent)
                 PS.<|>
                 fieldContent
                 


fieldContent :: GenParser Char st String
fieldContent = PS.many (noneOf "\t\n")


-- The end of line character is \n
eol :: GenParser Char st Char
eol = char '\n'
            
parseOM input = parse valouevFile "(unknown)" input



extract_frags :: (String, String, String, [Float]) -> [Float]
extract_frags (_, _, _, frags) = frags

frag_delim :: [Float]
frag_delim = [100000.0]

dumpnode :: (Float, Int) -> Put
dumpnode (label, value) = do
  putWord32host $ round label
  putWord32host $ fromIntegral value

enumerate :: [Float] -> [(Float, Int)]
enumerate order_list = (zip order_list [1..])
                       
dumpOrderList :: [Float] -> Put
dumpOrderList order_list = do
  let enumerated = enumerate order_list
  let actions = fmap dumpnode enumerated   -- broken, need to do full list, not just first element
  sequence actions
  return ()

         
dumpNodes :: [[Float]] -> Put
dumpNodes skip_nodes = do
  let num_nodes = fromIntegral $ sum $ fmap length skip_nodes
  putWord32host $ fromIntegral $ num_nodes + 2 -- for initial and final
  dumpnode initnode
  let actions = fmap dumpOrderList skip_nodes
  sequence actions
  dumpnode finalnode
  return ()

 -- this attaches numbers in the same order as they are written to disk, FIXME: figure out a way to assign these numbers at the time they are written 
enumerateOrderList :: (Int, [a]) -> [(a, Int)]
enumerateOrderList (start, list) = zip list [start..]

enumerateNodes :: [[a]] -> [[(a, Int)]]
enumerateNodes skip_nodes = let lengths = fmap length skip_nodes
                                prefixes = scanl (+) 1 lengths
                                prefix_list_pairs = zip prefixes skip_nodes in
                            fmap enumerateOrderList prefix_list_pairs
                                
max_skipnode = 3
initnode = (0.0, 2100000000) -- fixme, should be 2**32 - 1 but not sure if unsigned clean
finalnode = (0.0, 0)

--FIXME Edge and Node should probably be types
dumpEdgePiece :: (Float, Int) -> Put
dumpEdgePiece (label, pos) = do putWord32host $ fromIntegral pos

dumpEdge :: [(Float, Int)] -> Put
dumpEdge enum_nodes = do let actions = fmap dumpEdgePiece enum_nodes
                         sequence actions
                         return ()
                                                                   
dumpEdges :: [[(Float, Int)]] -> Put
dumpEdges edge_list = do putWord32host $ fromIntegral $ length edge_list
                         sequence $ fmap dumpEdge edge_list
                         return ()

dumpGraph all_skipnodes edges = do dumpNodes all_skipnodes
                                   dumpEdges edges -- FIXME: TODO

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
   Right x -> BL.hPut ohdl $ runPut $ dumpGraph all_skipnodes edges
              where  nodes = intercalate frag_delim $ fmap extract_frags x
                     skipnode_list n = nth_skipnodes n nodes
                     all_skipnodes = fmap skipnode_list [max_skipnode..1]
                     all_enumerated_skipnodes = enumerateNodes all_skipnodes
                     num_all_nodes = length $ concat all_enumerated_skipnodes
                     enumerated_initnode = (0.0, 0)
                     enumerated_finalnode = (0.0, num_all_nodes)
                     junction_nodes_pair = zip ([[enumerated_initnode]] ++  (rights all_enumerated_skipnodes)) ((lefts all_enumerated_skipnodes) ++ [[enumerated_finalnode]]) -- FIXME need start and end nodes
                     product (left, right) = sequence left right
                     products a = concat $ map product a 
                     edges = concat $ fmap product  junction_nodes_pair
  
  hClose ohdl  


-- can probably generate the middle portion of the left and right ends of edges using: take, drop, and transpose
-- 'sequence' will form the cartesian product of lists of lists

-- takes a list of skip lists and returns a list of heads; assumes shortest list first : FIXME: check for or eliminate this assumption
-- lefts means all the skipnodes have their left edge aligned, which actually means they form the right side of a junction
lefts :: [[a]] -> [[a]]
lefts [] = []
lefts [[]] = []
lefts ([]:skip_lists) = fmap head skip_lists : lefts ( fmap tail skip_lists )
lefts skip_lists = fmap head skip_lists : lefts ( fmap tail skip_lists )

-- FIXME: rights is inefficient, as it forces the construction and retension of intermediate lists, this could be done lazily and allow the used heads to be GC'd
rights :: [[a]] -> [[a]]
rights skip_lists = reverse $ lefts $ fmap reverse skip_lists
-- thanks to http://stackoverflow.com/questions/27726739/implementing-an-efficient-sliding-window-algorithm-in-haskell
windows :: Int -> [a] -> [[a]]
windows m = transpose . take m . tails

nth_skipnodes n nodes  = fmap sum $ windows n nodes

