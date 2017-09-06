awk '/Variation_project_Polymorphism/' WS254.gff2 > WS254.variants.gff2
awk '{ if ( $1 == "CHROMOSOME_II" ) { print } }' WS254.variants.gff2 | awk -F "\t" '{ if ( $4 <= 13098000 && $4 >= 12944000 ) { print } }' > WS254.CHII.NIL59.Variants.gff2
