library("tidyverse")

results = read.csv("results.csv")

resultspartialfiltered = filter(results, time_for_database < 20) #remove refhashtable

resultsfiltered = filter(resultspartialfiltered, Bucket_Number_with_pow_2 > 9) #PLOTBUG with BNWP <= 9

selectedresults = select(resultsfiltered, Bucket_Number_with_pow_2, Genome0) #Just the genome0


ggplot(selectedresults, aes(Bucket_Number_with_pow_2, Genome0, group = 1)) +
  geom_line(color="navyblue") +
  geom_hline(yintercept=0.228, size=0.25, color="seagreen") +
  coord_cartesian(xlim=c(9,26)) +
  theme_linedraw(base_line_size = 0.5) +
  labs(
        x="Bits (for 2 ^ bits Buckets number)",
        y="Jaccard Index (of the query correponding genome)",
        title="Precision of index table in function of buckets number",
        subtitle=" with the reference of the hashtable",
        caption="Experimental condition : 1 query (2.3Mo) among 3 genomes (30Mo)") +
  annotate (geom="text",x = 11, y = 0.238, label = "Absolute Jaccard Index", color="seagreen")


ggsave("plot.png", width = 7, height = 7)
