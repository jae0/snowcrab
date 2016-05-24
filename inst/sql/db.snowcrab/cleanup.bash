

#  scs_main
gunzip scs_main.gz
gawk 'BEGIN {RS="\n"; FS=";"; OFS=";"} \
      NF>5 && NR>2 && !/LOADING/ && !/----/ && !/SQL/ && !/select/ \
      {gsub(/[ ]/,""); print $0}'  scs_main > scs_main.out
gzip scs_main


#  scs_area
gunzip scs_area.gz
gawk 'BEGIN {RS="\n"; FS=";"; OFS=";"} \
      NF>5 && NR>2 && !/TRAWL/ && !/----/ && !/SQL/ && !/select/ \
      {gsub(/[ ]/,""); print $0}'  scs_area > scs_area.out
gzip scs_area


#  scs_count ... 
gunzip scs_count.gz
gawk 'BEGIN {ORS="~"} !/LOADING/ && !/TMM/ && !/----/ && !/SQL/ && !/select/ \
      {gsub(/[ ]/,""); print $0}'  scs_count > tmp
# use sed to change record separators to newlines and remove the first ";"
sed -e 's/~~/\n/g' -e 's/~/;/g' -e 's/^;//1'  tmp > scs_count.out
gzip scs_count



