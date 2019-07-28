library(sets)
path_to_data <- "competitor_pack/data/fire-data.csv"
data = read.csv(path_to_data, nrows = 100)
q=set("1", "10", "10,11", "10,14", "11", "11,13", "11,14", "11,15", "12","13", "13,5", "13,6", "14", "14,3", "14,4", "15", "16", "17", "18","19", "2", "2,5", "2,6", "20", "21", "22", "23", "24", "25", "26","27", "28", "29", "29,3", "3", "3,31", "30", "31", "32", "4", "4,7","4,9", "5", "5,6", "6", "6,7", "7", "8", "9")

print(names(data))
queries<-c(c())#
for (marginal in q){
    print(marginal)
    marg=as.set(strsplit(marginal,',')[[1]])
    m_list <- c()#
    for (attr in marg){
     #   print(names(data)[as.numeric(attr)])
        m_list = c(m_list,paste(names(data)[as.numeric(attr)]))    #  
    }
    print(length(m_list))
    queries=c(queries,paste(c(m_list)))#
}
print("ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss")
print(queries)
print(length(queries))
