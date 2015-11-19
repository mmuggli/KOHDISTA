import Text.ParserCombinators.Parsec

valouevFile :: GenParser Char st [(String, String, String, [Float])]
valouevFile =
            do  result <- many record
                eof
                return result

record :: GenParser Char st (String, String, String, [Float])
record = do map_name <- fieldContent
            eol
            char '\t'
            enz_name <- fieldContent
            char '\t'
            enz_acr_name <- fieldContent
            char '\t'
            frags <- sepBy fieldContent (char '\t')
            eol
            eol
            return (map_name, enz_name, enz_acr_name, fmap floatify frags)




fieldContent :: GenParser Char st String
fieldContent = many (noneOf "\t\n")


-- The end of line character is \n
eol :: GenParser Char st Char
eol = char '\n'
            
parseOM input = parse valouevFile "(unknown)" input

floatify :: String -> Float
floatify s = read s

