############################################################################
#################  METADATA   ##############################################
############################################################################

# mega.phy: time-calibrated molecular phylogeny of 13 000+ plant species 
            # (Zane et al. 2014. Nature; 
            # updated by Qian and Jin. 2015. Journal of PLant Ecology)
# comm1: pseudo community matrix that contains spp of interest
# sample: list of binomial latin names of all spp of interest 
# phylo: pruned subtree that contains spp of interest

############################################################################

#### Load Libraries ####
library(ape)
library(picante)

#### Load megaphylogeny.new ####

mega.phy <- read.tree(file.choose())

### Create pruned phylogeny from species of interest ###

# Species of interest
# NOTE: Rhodiola_rhodantha acts as a congeneric relative (proxy) for Rhodiola_rosea
sample <- c("Empetrum_nigrum",
            "Gaultheria_procumbens",
            "Vaccinium_vitis-idaea",
            "Danthonia_spicata",
            "Deschampsia_flexuosa",
            "Poa_compressa",
            "Campanula_rotundifolia",
            "Plantago_maritima",
            "Sedum_acre",
            "Solidago_bicolor",
            "Sagina_procumbens",
            "Rhodiola_rhodantha",
            "Sedum_spurium")

# Create a pseudo community matrix
comm1 <- rbind(sample)
colnames(comm1) <- comm1

# Prune the megaphylogeny 
phylo <- prune.sample(comm1, mega.phy)

# Rename Rhodiola_rhodantha as "Sedum_rosea" in the phylogeny
# Double-check indexing; 5 could be some other species 
phylo$tip.label[5] <- "Sedum_rosea"


############################################################################

# FINAL RESULT: Plot the pruned phylogeny
plot.phylo(phylo, label.offset = 2)
axisPhylo()
