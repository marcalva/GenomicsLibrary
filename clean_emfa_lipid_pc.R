#
# Handle EMFA and lipid data

emfa_file_name = "work/EMFA/metsim793_ffa_ery_only.csv"
alltraits_file_name = "work/Mete0914/METSIM-Baseline_and_FollowUp_FatBiopsySamples_07012014.csv"
# metsim494_ids_fn = "work/metsim/metsim494_vb_unique_match.txt"

emfa = read.csv(emfa_file_name, header=TRUE, row.names=1, stringsAsFactors=F, check.names=F)
alltraits = read.csv(alltraits_file_name, header=TRUE, row.names=1, stringsAsFactors=F, check.names=F, na.string="#NULL!")
# samples = read.table(metsim494_ids_fn, header=TRUE, stringsAsFactors=F, check.names=F)
# samples = as.character(samples[,2])
samples = rownames(emfa)[rownames(emfa) %in% rownames(alltraits)]
emfa = emfa[samples,]
alltraits = alltraits[samples,]


###################################
# EMFA
###################################

emfa_er = as.matrix(emfa[,c(seq(from=1,to=46,by=2), seq(from=49,to=58,by=2))])
# emfa_er_norm = apply(emfa_er, 2, function(x) qqnorm(resid(summary(lm(x ~ alltraits$Age, na.action=na.exclude))), plot.it=FALSE)$x)
emfa_er_norm = apply(emfa_er, 2, function(x) resid(summary(lm( log(x + .5) ~ alltraits$Age , na.action=na.exclude))) )

###################################
# Metabolic Traits
###################################
mtNames = c("S_hdlc", "S_ldlc", "S_tottg", "S_totalc", "BMI")
mt = alltraits[,mtNames]
# mt_norm = apply(mt , 2, function(x) qqnorm(resid(summary(lm(x ~ alltraits$Age, na.action=na.exclude))), plot.it=FALSE)$x)
mt_norm = apply(mt , 2, function(x) resid(summary(lm( log(x + .5) ~ alltraits$Age , na.action=na.exclude))) )


###################################
# Amino Acids
###################################
aaNames = c("Gln", "Gly", "His", "Ile", "Leu", "Phe", "Tyr", "Val")
aa = alltraits[,aaNames]
# aa_norm = apply(aa, 2, function(x) qqnorm(resid(summary(lm(x ~ alltraits$Age, na.action=na.exclude))), plot.it=FALSE)$x)
aa_norm = apply(aa, 2, function(x) resid(summary(lm( log(x + .5) ~ alltraits$Age , na.action=na.exclude)) ))


###################################
# Combine
###################################


emfaAaMt = cbind(emfa_er_norm, mt_norm, aa_norm)

# Remove NAs
emfaAaMtNoNA = emfaAaMt[rowSums(is.na(emfaAaMt)) == 0,]

emfaAaMtNoNAPc = prcomp(t(emfaAaMtNoNA), retx=TRUE, scale = TRUE)

# Create factors for cluster
traitCategories = c(rep("Saturated", 7), rep("Mono_Unsaturated", 6), rep("Poly_Unsaturated", 10), rep("Ratio", 5), rep("Metabolic", 5), rep("Amino_Acid", 8))
traits = factor(traitCategories)

library(ggplot2)
ggplot(as.data.frame(emfaAaMtNoNAPc$x), aes(x=PC1, y=PC2, colour=traits)) + geom_point()

emfaAa = cbind(emfa_er_norm, aa_norm)
emfaAaNoNA = emfaAa[rowSums(is.na(emfaAa)) == 0,]

emfaAaNoNAPc = prcomp(t(emfaAaNoNA), retx=TRUE, scale = TRUE)
traitCategories = c(rep("Saturated", 7), rep("Mono_Unsaturated", 6), rep("Poly_Unsaturated", 10), rep("Ratio", 5), rep("Amino_Acid", 8))
traits = factor(traitCategories)
ggplot(as.data.frame(emfaAaNoNAPc$x), aes(x=PC1, y=PC2, colour=traits)) + geom_point()



