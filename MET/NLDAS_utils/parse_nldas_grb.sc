#mkdir -p MET_month_files
for i in *.grb; do
    [ -f "$i" ] || break
    wgrib $i -v | wgrib -text -i $i -o $i.txt.temp > metadata.txt
    awk '$0!="1 1" {print}' $i.txt.temp > $i.txt
    #mv $i.txt ./MET_month_files 
    rm $i.txt.temp
    #mv metadata.txt ./MET_month_files
done

wgrib -v $i > metadata_full.txt
#cp metadata_full.txt ./MET_month_files
