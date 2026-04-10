library(ggplot2)

args = commandArgs(trailingOnly=TRUE)

if (length(args)!=2) {
  stop("Usage: Rscript --vanilla run_plot_metagene_density.R Input(metagene.txt) Output(metagene.pdf)", call.=FALSE)
}

dat <- read.table(args[1],header=T,sep="\t")

d <- ggplot(data=dat, aes(x=Rel_Location, group=Sample, color=Sample)) + geom_density(adjust=1.0)+theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + xlab("5UTR_CDS_3UTR") + ylab("Peak density")

ggsave(plot=d,width=6,height=3,dpi=300,filename=args[2])
