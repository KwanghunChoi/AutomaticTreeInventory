# This is the code to create inventory automatically with pre-processed images from YOLOv3 and semantic segmentation with cityscape pretrained model.

AutomaticCode <-
    function(RawImage, YoloImage, YoloTable, SemsegImage, processDir, VH=2.5){
        # Install and load essential packages for analysis
        pacman::p_load("imager", "tidyverse", "foreach", "segmented", "rmarkdown", "pander", "strucchange")

        # Extract Metadata from the file name of the Raw Image
        FileName <- gsub(".jpg", "", tail(unlist(strsplit(RawImage, "/")), 1))
        MetaData <- data.frame(t(unlist(strsplit(FileName, "_"))))
        colnames(MetaData) <- c("Location", "Year", "ID", "Latitude", "Longitude", "Heading")

        # Load images
        Orig_img <- load.image(RawImage)     # original image for image size
        Yolo_img <- load.image(YoloImage)     # original image for image size
        SmSg_img <- load.image(SemsegImage)   # Image from semantic segmentation

        # Load location of trees
        # (Bounding box from YOLO (left top right bottom class))
        Yolo_bbox <- read.csv(YoloTable)
        # Information of original pictures
        W <- nrow(Orig_img) # Width of the picture
        H <- ncol(Orig_img) # Height of the picture
        VP <- (H+1)/2       # View point of the picture (!H+1 is required to point out the middel line of the image)
        
        # Required functions
        # Split image with imsub
        HorSpliter <- function(Img, BBox){ Img %>% imsub(x >= BBox[1] & x <= BBox[3])}
        VerSpliter <- function(Img, BBox){ Img %>% imsub(y >= BBox[2] & y <= BBox[4])}

        # Image extractor
        ## Supl. Extract colors of vegetation (107, 142, 35) from panoptic-deeplab hompage (https://github.com/bowenc0221/panoptic-deeplab/issues/42)
        TreeExtractor <- 
            function(Img){
                IO <- Img
                R(IO)[which(!R(IO)[] * 255 == 107)] <- NA
                G(IO)[which(!G(IO)[] * 255 == 142)] <- NA
                B(IO)[which(!B(IO)[] * 255 == 35)] <- NA
                return(IO)
            }

        GrndExtractor <- 
            function(Img){
                Road <- Img
                R(Road)[which(!R(Road)[] * 255 == 128)] <- 0
                G(Road)[which(!G(Road)[] * 255 == 64)] <- 0
                B(Road)[which(!B(Road)[] * 255 == 128)] <- 0
                SideWalk <- Img
                R(SideWalk)[which(!R(SideWalk)[] * 255 == 244)] <- 0
                G(SideWalk)[which(!G(SideWalk)[] * 255 == 35)] <- 0
                B(SideWalk)[which(!B(SideWalk)[] * 255 == 232)] <- 0
                Terrain <- Img
                R(Terrain)[which(!R(Terrain)[] * 255 == 152)] <- 0
                G(Terrain)[which(!G(Terrain)[] * 255 == 251)] <- 0
                B(Terrain)[which(!B(Terrain)[] * 255 == 152)] <- 0
                IO <- Road + SideWalk + Terrain
                R(IO)[which(R(IO)[] * 255 == 0)] <- NA
                G(IO)[which(G(IO)[] * 255 == 0)] <- NA
                B(IO)[which(B(IO)[] * 255 == 0)] <- NA
                return(IO)
            }

        # Detecting base point of tree
        BaseDetector <- 
            function(Img, Grnd){
                # Extract the range of tree height 
                TreeOnly_col <- as.numeric(colSums(R(Img), na.rm=T)) %>% ifelse(. < 1, NA, .)
                TreeOnly_min <- ((TreeOnly_col * 0 + 1) * c(1:H)) %>% min(.,na.rm=T)
                TreeOnly_max <- ((TreeOnly_col * 0 + 1) * c(1:H)) %>% max(.,na.rm=T)
                GrndOnly_col <- as.numeric(colSums(R(Grnd), na.rm=T)) %>% ifelse(. < 1, NA, .)
                GrndOnly_min <- ((GrndOnly_col * 0 + 1) * c(1:H)) %>% min(.,na.rm=T)
                Base         <- max(TreeOnly_max, GrndOnly_min)
                Diff         <- Base - TreeOnly_max + 1         # Plus one is required for including base cells
                return(c(Base, Diff, TreeOnly_min, TreeOnly_max))
            }

        # Image analysis tools
        StructureAnalyzer <- 
            function(out, H){
                # Remove the row with zero vegetation pixels
                out.nonzero <- out[out[]>0]
                out.df      <- data.frame(PixW = out.nonzero)
                out.df$PixH <- c(1:nrow(out.df)) 
                cellNumber  <- as.numeric(gsub("X", "", row.names(out.df)))

                # Check the first change point of the plot
                out.cp <- strucchange::Fstats(PixW ~ 1, data=out.df)

                ## 3. Create Tree inventory table
                # Create Tree inventory in pixels
                VP       <- (H + 1) / 2
                Tree_W   <- max(out.df$PixW, na.rm=T)                                           # Crown width
                CellToVP <- ceiling(abs(cellNumber - VP))
                Tree_H   <- sum(1/cos(asin(CellToVP/(H/2))), na.rm=T)                           # Tree Height
                BreakPt  <- out.cp$breakpoint
                Tree_BH  <- sum(1/cos(asin(CellToVP[1:BreakPt]/(H/2))), na.rm=T)                # Height below crown.
                # DBH which is assumed to be a median width below crown.
                Tree_DBH <- quantile(out.df$PixW[1:Tree_BH], 0.1) 

                # Create Table for the metric.
                M <- data.frame(Tree_H, Tree_W, Tree_BH, Tree_DBH, TreeTop = min(cellNumber))
                return(M)
            }

        # Combined function
        # VH: height of the camera
        TreeInventory  <- 
            function(Img, VH = 2.5, LocInfo = TreeLoc){ # VH of GSV is known as 2.5m 8.2 feet
                W <- nrow(Img) # Width of the picture
                H <- ncol(Img) # Height of the picture
                VP <- (H+1)/2      # View point of the picture
                BM <- foreach(i = 1:nrow(LocInfo), .combine=rbind) %do% {
                    tryCatch({
                        HorSplitImg     <- HorSpliter(Img, as.numeric(LocInfo[i,]))
                        IndTreeOnly     <- TreeExtractor(HorSplitImg)
                        Grnd            <- GrndExtractor(HorSplitImg)
                        BaseDetector    <- BaseDetector(IndTreeOnly, Grnd)
                        Base            <- BaseDetector[1]
                        Diff            <- BaseDetector[2]              # I can calculate the corrected diff but it might be almost same.
                        TreeTop         <- BaseDetector[3]              # I don't know that I really need this.
                        TreeBottom      <- BaseDetector[4]
                        # Calculate the pixel width of the tree
                        BaseToVP        <- ceiling(abs((Base - VP)))    # As they have less distortion we put them just like that. (Basic Pixel width)
                        CorrectedPixSum <- sum(1/cos(asin(c(1:BaseToVP)/(H/2))), na.rm=T)
                        PixelWidth      <- VH/CorrectedPixSum
                        TreeFit         <- VerSpliter(IndTreeOnly, as.numeric(LocInfo[i,]))
                        PixelCnt        <- rev(sapply(data.frame(ceiling(TreeFit[,,1])), FUN = function(x){sum(x, na.rm=T)}))
                        names(PixelCnt) <- paste0("X", c(LocInfo$bottom[i]:LocInfo$top[i]))
                        PixelStr        <- StructureAnalyzer(PixelCnt, H)
                        CorrectedDiff   <- ifelse(Diff == 0, 0, sum(1/cos(asin(ceiling(abs((c(TreeBottom:Base) - VP)))/(H/2))), na.rm=T))
                        result          <- (PixelStr %>% mutate(Tree_H = Tree_H + CorrectedDiff, Tree_BH = Tree_BH + CorrectedDiff)) * PixelWidth 
                        result$TreeTop  <- PixelStr$TreeTop
                        result$ID       <- i
                        return(result)
                    }, error=function(e){})
                    }
            }

        # structure of the tree
        Tree.df <- TreeInventory(Img = SmSg_img, VH = 2.5, LocInfo = Yolo_bbox) 

        # Location of the tree
        Tree.Inv <- data.frame(Tree.df, Yolo_bbox[Tree.df$ID,])

        Tree.Inv <- Tree.Inv %>% 
                    # The field of view (FOV) of Google street view is known as 127 degree and 63.5 is the half of the FOV.
                    mutate(Dist  = abs(Tree_H - VH)/tan(asin(ceiling(abs(TreeTop-(H+1)/2))/(H/2)))) %>%
                    # Horizontally, 1 pixel denotes the W/360 degree
                    mutate(Angle = (as.numeric(MetaData$Heading) + (0.5 * (left + right) - W/2) * 360/W) %%360) %>%
                    # Polar coordinate to cartesian coordinate
                    mutate(rel_x = Dist * cos(Angle), rel_y = Dist * sin(Angle)) %>%
                    mutate(c_lat = as.numeric(MetaData$Latitude), c_lon = as.numeric(MetaData$Longitude))
        write.csv(Tree.Inv, paste0(processDir, FileName, "_result.csv") , row.names=F)
        return(Tree.Inv)
    }
