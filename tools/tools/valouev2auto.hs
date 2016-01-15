import Debug.Trace
import Text.ParserCombinators.Parsec as PS
import System.IO
import System.Environment
import Data.List
import qualified Data.ByteString.Lazy as BL
import Data.Binary.Put
--import Control.Applicative as A
bin_size :: Int
bin_size = 20

-- the order of the automaton.  This counts the backbone as a degenerate case of skip nodes, you might call it the 0th order skip nodes.  So an automaton with a backbone and skipnodes that sum two consecutive nodes in the backbone for each skipnode would have a value of 2.                    
max_skipnode :: Int                                
max_skipnode = 3
               
desorption_thresh :: Float
desorption_thresh = 1.0
                    
valouevFile :: GenParser Char st [(String, String, String, [Float])]
valouevFile =
            do  result <- PS.many record
                eof
                return result

record :: GenParser Char st (String, String, String, [Float])
record = do map_name <- fieldContent ; _ <- eol
            enz_name <- parseEnzyme ; _ <- char '\t'; enz_acr_name <- fieldContent; _ <- char '\t'; frags <- sepBy fieldContent (char '\t'); _ <- eol
            _ <- eol
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

parseOM :: [Char] -> Either ParseError [(String, String, String, [Float])]  
parseOM input = parse valouevFile "(unknown)" input



extract_frags :: (String, String, String, [Float]) -> [Float]
extract_frags (_, _, _, frags) = frags

frag_delim :: [Float]
frag_delim = [100000.0]

dumpnode :: (Int, Int) -> Put
dumpnode (lab, value) = do if lab == 0
                           then trace (show value) $ putWord32host $ fromIntegral lab
                           else  putWord32host $ fromIntegral lab
                           putWord32host $ fromIntegral value

enumerate :: [a] -> [(a, Int)]
enumerate list = (zip list [1..])
                       
dumpOrderList :: [Int] -> Put
dumpOrderList order_list = do
  let enumerated = enumerate order_list
  let actions = fmap dumpnode enumerated   -- broken, need to do full list, not just first element
  sequence_ actions

make_skipnode n = if even n
                  then n + 1
                  else n
                       
make_backbone n = if odd n
                  then n - 1
                  else n

process_backbone_node :: Float -> Int                  
process_backbone_node lab = if retval == 0
                        then trace ("proc_backbone: " ++ (show lab) )  retval
                        else retval
                            where retval :: Int
                                  retval =  fromIntegral $ make_backbone $  round (lab * 1000)
                        

process_backbone_nodes :: [Float] -> [Int]
process_backbone_nodes labs = fmap process_backbone_node labs

process_actual_skipnode :: Float -> Int
process_actual_skipnode = (\lab -> fromIntegral $ make_skipnode $ round (lab * 1000))                              
process_actual_skipnodes :: [Float] -> [Int]                             
process_actual_skipnodes labs = fmap process_actual_skipnode labs
         
--dumpNodes :: [[Float]] -> Put
dumpNodes skip_nodes nodes = do
  let num_nodes = fromIntegral $ sum $ fmap length skip_nodes
  putWord32host $ fromIntegral $ num_nodes + 2 -- for initial and final
  --initnode :: (Int, Int)
  let initnode = (2100000000, 0) -- fixme, should be 2**32 - 1 but not sure if unsigned clean

  --finalnode :: (Int, Int)
  let finalnode = (process_backbone_node 0.0, (length nodes) + 1)

  dumpnode initnode
  let actual_skipnodes = init  skip_nodes
  let backbone = [last skip_nodes]
  let actions = (fmap dumpOrderList (fmap process_actual_skipnodes actual_skipnodes)) ++ (fmap dumpOrderList (fmap process_backbone_nodes backbone))
  sequence_ actions
  dumpnode finalnode
  --return ()

 -- this attaches numbers in the same order as they are written to disk, FIXME: figure out a way to assign these numbers at the time they are written 
enumerateOrderList :: (Int, [a]) -> [(a, Int)]
enumerateOrderList (start, list) = zip list [start..]

enumerateNodes :: [[a]] -> [[(a, Int)]]
enumerateNodes skip_nodes = let lengths = fmap length skip_nodes
                                prefixes = scanl (+) 1 lengths
                                prefix_list_pairs = zip prefixes skip_nodes in
                            fmap enumerateOrderList prefix_list_pairs


quantize :: Int -> Int
quantize val = if (val `mod` bin_size) < (bin_size `quot` 2)
               then val - (val `mod` bin_size)
               else val - (val `mod` bin_size) + bin_size

float_quantize :: Float -> Float
float_quantize val = (fromIntegral (quantize (round (1000.0 * val)))) / 1000.0


--FIXME Edge and Node should probably be types
dumpEdgePiece :: (Float, Int) -> Put
dumpEdgePiece (_, pos) = do putWord32host $ fromIntegral pos

dumpEdge :: [(Float, Int)] -> Put
dumpEdge enum_nodes = do let actions = fmap dumpEdgePiece enum_nodes
                         sequence_ actions
                         --return ()
                                                                   
dumpEdges :: [[(Float, Int)]] -> Put
dumpEdges edge_list = do putWord32host $ fromIntegral $ length edge_list
                         sequence_ $ fmap dumpEdge edge_list
                         --return ()

--dumpGraph :: [[Float]] -> [[(Float, Int)]] -> Put
dumpGraph all_skipnodes edges nodes = do dumpNodes all_skipnodes nodes
                                         dumpEdges edges -- FIXME: TODO


-- FIXME: enable this function
add_reversed :: [[a]] -> [[a]]
add_reversed [] = []
add_reversed a = (head a) : (reverse (head a)) : (add_reversed (tail a))

-- FIXME: trim tips of rmaps

main :: IO ()
main = do
  args <- getArgs
  let fname = (head args)
  let ofname = (last args)
  hdl <- openFile fname ReadMode 
  ohdl <- openFile ofname WriteMode
  contents <- hGetContents hdl
  -- case parseOM contents of
  --  Left x -> print $ show $ x
  --  Right x -> print $ show $ intercalate frag_delim $ fmap extract_frags x
  case parseOM contents of
   Left x -> print $ show $ x
   Right x -> sequence_ [show_stats,  dump_file ]
              where  show_stats = print $ "nodes: " ++ (show num_all_nodes) ++ " edges: " ++ (show (length edges)) ++ " skip edge bundles: " ++ (show (length skip_edge_list))
                     dump_file = BL.hPut ohdl $ runPut $ dumpGraph all_skipnodes edges nodes
                     nodes = {-filter (\a -> a > 0.001) $-} fmap float_quantize $ intercalate frag_delim $ {- add_reversed $ -}  fmap extract_frags x
                     skipnode_list n = take ((length nodes) - (n - 1)) $ nth_skipnodes n nodes
                     all_skipnodes = fmap skipnode_list $ reverse [1..max_skipnode]
                     all_enumerated_skipnodes = enumerateNodes all_skipnodes -- this enumeration is file position, not backbone position
                     num_all_nodes = length $ concat all_enumerated_skipnodes
                     enumerated_initnode = (0.0, 0)
                     enumerated_finalnode = (0.0, num_all_nodes + 1)
                     
                     all_rights :: [[(Float, Int)]]
                     all_rights = [[enumerated_initnode]] ++ (rights all_enumerated_skipnodes)

                     all_lefts :: [[(Float, Int)]]
                     all_lefts = (lefts all_enumerated_skipnodes) ++ [[enumerated_finalnode]]

                     junction_nodes_pair :: [     ([(Float, Int)], [(Float, Int)])    ]
                     junction_nodes_pair = zip all_rights all_lefts

                     product ::  ([(Float, Int)], [(Float, Int)])   -> [[(Float, Int)]]
                     product (left, right) =  sequence [left, right]

                     --products a = concat $ map product a 
                     edge_lists :: [[[(Float, Int)]]]
                     edge_lists = fmap product  junction_nodes_pair
                     edge_skipable :: [([(Float, Int)], [(Float, Int)])]
                     edge_skipable = zip (all_rights ++ [[enumerated_finalnode]]) ([[enumerated_initnode]] ++ all_lefts)
                     actual_skip = skip_edge_pairs edge_skipable
                     skip_edge_list = fmap product $ fmap productable actual_skip
                     edges :: [[(Float, Int)]]
                     edges =  (concat  edge_lists) ++ (concat skip_edge_list) 
  
  hClose ohdl  

productable :: (([(Float, Int)], [(Float, Int)]), ([(Float, Int)], [(Float, Int)])) -> ([(Float, Int)], [(Float, Int)])
productable ((a,b), (c,d)) = (a,d)

skip_edge_pairs :: [([(Float, Int)], [(Float, Int)])] -> [(([(Float, Int)], [(Float, Int)]), ([(Float, Int)], [(Float, Int)]))]
skip_edge_pairs (x:y:[]) = []
skip_edge_pairs (x:xs) = if small $ head xs
                         then skip_edge_run x xs
                         else skip_edge_pairs xs

skip_edge_run :: ([(Float, Int)], [(Float, Int)]) -> [([(Float, Int)], [(Float, Int)])] -> [(([(Float, Int)], [(Float, Int)]), ([(Float, Int)], [(Float, Int)]))]
skip_edge_run x (y:[]) = []
skip_edge_run x xs = if small $ head xs
                     then (x, head (tail xs)) : skip_edge_run x (tail xs)
                     else skip_edge_pairs xs

small :: ([(Float, Int)], [(Float, Int)]) -> Bool
small (rights, lefts) =  case (last rights) of (label, pos) ->  label < desorption_thresh

-- can probably generate the middle portion of the left and right ends of edges using: take, drop, and transpose
-- 'sequence' will form the cartesian product of lists of lists

-- takes a list of skip lists and returns a list of heads; assumes shortest list first and decreasing size : FIXME: check for or eliminate this assumption
-- lefts means all the skipnodes have their left edge aligned, which actually means they form the right side of a junction
lefts :: [[a]] -> [[a]]
lefts [] = []
lefts [[]] = []
lefts ([]:skip_lists) =  fmap head skip_lists : lefts ( fmap tail skip_lists )
lefts skip_lists =  fmap head skip_lists : lefts ( fmap tail skip_lists )

-- FIXME: rights is inefficient, as it forces the construction and retension of intermediate lists, this could be done lazily and allow the used heads to be GC'd
rights :: [[a]] -> [[a]]
rights skip_lists = reverse $ lefts $ fmap reverse skip_lists
-- thanks to http://stackoverflow.com/questions/27726739/implementing-an-efficient-sliding-window-algorithm-in-haskell
windows :: Int -> [a] -> [[a]]
windows m = transpose . take m . tails

nth_skipnodes :: Num b => Int -> [b] -> [b]
nth_skipnodes n nodes  = fmap sum $ windows n nodes

--TODO: mark backbone specially
