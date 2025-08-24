# Clear environment and set seed for reproducibility 
rm(list = ls()) 
set.seed(33270961)

# ============================================================================== 
# LOAD REQUIRED LIBRARIES 
# ============================================================================== 

library(tm)                # Text mining 
library(SnowballC)         # Stemming 
library(proxy)             # Distance calculations 
library(SentimentAnalysis) # Sentiment analysis 
library(igraph)            # Network analysis 
library(RColorBrewer)      # Color palettes

# ============================================================================== 
# DATA LOADING AND PREPROCESSING 
# ============================================================================== 

# Set path to corpus folder 
cname <- file.path(".", "corpus") 
# Read documents into a Corpus object 
docs <- Corpus(DirSource(cname)) 
print(summary(docs)) 
# Define common words to remove (in addition to standard stopwords) 
common_words <- c( 
  "also", "back", "can", "come", "even", "first", "get", "just", "last", "like", 
  "new", "now", "said", "take", "time", "will", "make", "want", "year", "many", 
  "much", "thing", "way", "one", "two", "three", "should", "would", "may", "might",  
  "say", "around", "still", "call", "look", "need", "made", "before", "because",  
  "include", "since", "see", "move", "meet", "give", "only", "think", "become", 
  "know", "though", "day", "week", "very", "yet", "part", "tell", "any", "almost", 
  "always", "set", "turn", "happen", "recent", "far", "past", "return", "follow", 
  "during", "start", "end", "try", "clear" 
) 
# Create regex pattern for common word removal 
pattern <- paste0("\\b(", paste(common_words, collapse = "|"), ")\\b") 
# Text preprocessing pipeline 
# Remove numbers 
docs <- tm_map(docs, removeNumbers) 
# Remove punctuation 
docs <- tm_map(docs, removePunctuation)  
# Remove quotes 
docs <- tm_map(docs, content_transformer(function(x) gsub("['‘’“”\"]+", " ", x))) 
# Convert to lowercase 
docs <- tm_map(docs, content_transformer(tolower)) 
# Remove common words 
docs <- tm_map(docs, content_transformer(function(x) gsub(pattern, " ", x))) 
# Remove English stopwords 
docs <- tm_map(docs, removeWords, stopwords("english")) 
# Remove extra whitespace 
docs <- tm_map(docs, stripWhitespace) 
# Stem words 
docs <- tm_map(docs, stemDocument, language = "english")

# ============================================================================== 
# DOCUMENT-TERM MATRIX CREATION 
# ============================================================================== 
                                    
# Create the Document-Term Matrix 
dtm <- DocumentTermMatrix(docs) 
# Remove sparse terms (keeping terms that appear in at least 30% of documents) 
dtm_sparse <- removeSparseTerms(dtm, 0.7) 
print(dim(dtm_sparse)) 
print(inspect(dtm_sparse)) 
# Convert to matrix and select 30 most informative tokens 
dtm_mat <- as.matrix(dtm_sparse) 
selected_tokens <- c( 
  "action", "countri", "former", "kill", "lead", "media", "play", "power", 
  "social", "state", "support", "war", "world", "decis", "forc", "issu", 
  "friend", "futur", "win", "announc", "defeat", "fight", "champion", 
  "final", "live", "fan", "event", "love", "great", "stori" 
) 
dtm_mat <- dtm_mat[, colnames(dtm_mat) %in% selected_tokens] 
View(dtm_mat)

# ============================================================================== 
# HIERARCHICAL CLUSTERING ANALYSIS 
# ==============================================================================
                                         
# Calculate cosine distance and perform hierarchical clustering 
cosine_dist <- proxy::dist(dtm_mat, method = "cosine") 
hc <- hclust(cosine_dist, method = "ward.D") 
# Plot dendrogram with 3 clusters 
plot(hc, hang = -1, main = "Dendrogram of Corpus (k = 3)") 
rect.hclust(hc, k = 3, border = "red") 
# Extract cluster assignments 
cluster_assignments <- cutree(hc, k = 3) 
print(cluster_assignments) 
# Create genre labels based on filename patterns 
namelist <- as.data.frame(as.table(summary(docs))) 
namelist <- head(namelist, length(docs)) 
genres <- ifelse(grepl("politics", namelist$Var1), "politics", 
                 ifelse(grepl("review", namelist$Var1), "review", "sport")) 
# Cross-tabulate genres vs clusters 
table(Genres = genres, Clusters = cluster_assignments)

# ============================================================================== 
# SENTIMENT ANALYSIS BY GENRE 
# ==============================================================================
                                         
# Perform sentiment analysis 
SentimentA <- analyzeSentiment(docs) 
SentimentA$Genre <- genres 
# Create side-by-side boxplots 
par(mfrow = c(1, 2)) 
boxplot(SentimentQDAP ~ Genre, data = SentimentA, frame = TRUE,  
        main = "SentimentQDAP by Genre") 
boxplot(PositivityQDAP ~ Genre, data = SentimentA, frame = TRUE,  
        main = "PositivityQDAP by Genre") 
# Display boxplot statistics 
boxplot(SentimentQDAP ~ Genre, data = SentimentA, plot = FALSE) 
boxplot(PositivityQDAP ~ Genre, data = SentimentA, plot = FALSE)                       
# Pairwise t-tests for positivity differences between genres 
# Test 1: Review vs Politics 
lhs <- SentimentA[SentimentA$Genre == "review", "PositivityQDAP"] 
rhs <- SentimentA[SentimentA$Genre == "politics", "PositivityQDAP"] 
t.test(lhs, rhs, alternative = "greater") 
# Test 2: Review vs Sport 
lhs <- SentimentA[SentimentA$Genre == "review", "PositivityQDAP"] 
rhs <- SentimentA[SentimentA$Genre == "sport", "PositivityQDAP"] 
t.test(lhs, rhs, alternative = "greater") 
# Test 3: Sport vs Politics 
lhs <- SentimentA[SentimentA$Genre == "sport", "PositivityQDAP"] 
rhs <- SentimentA[SentimentA$Genre == "politics", "PositivityQDAP"] 
t.test(lhs, rhs, alternative = "greater")

# ============================================================================== 
# DOCUMENT SIMILARITY NETWORK (SINGLE-MODE) 
# ============================================================================== 
                                         
# Convert DTM to binary (presence/absence of terms) 
dtm_bin <- (dtm_mat > 0) + 0 
# Compute shared term matrix (documents × documents) 
adj_mat <- dtm_bin %*% t(dtm_bin) 
# Remove self-loops (documents don't connect to themselves) 
diag(adj_mat) <- 0 
# Create igraph object for document similarity network 
g <- graph_from_adjacency_matrix(adj_mat, mode = "undirected", weighted = TRUE) 
# Reset plotting layout 
par(mfrow = c(1, 1)) 
# Plot basic document similarity network 
set.seed(33270961) 
plot(g, vertex.label.cex = 0.7, main = "Single-Mode Network: Shared Terms") 
# Calculate centrality measures (using inverse weights for distance-based measures) 
inv_weights <- 1 / E(g)$weight 
# Degree centrality: Number of direct connections per document 
deg <- degree(g) 
cat("Top degree centrality documents:\n") 
print(sort(deg, decreasing = TRUE)) 
# Betweenness centrality: Documents that bridge between groups 
betw <- betweenness(g, weights = inv_weights) 
cat("\nTop betweenness centrality documents:\n") 
print(sort(betw, decreasing = TRUE)) 
# Closeness centrality: How close each document is to all others 
close <- closeness(g, weights = inv_weights) 
cat("\nTop closeness centrality documents:\n") 
print(sort(close, decreasing = TRUE)) 
# Eigenvector centrality: Documents connected to other well-connected documents 
eig <- eigen_centrality(g)$vector 
cat("\nTop eigenvector centrality documents:\n") 
print(sort(eig, decreasing = TRUE)) 
# Community detection using Louvain method 
comm <- cluster_louvain(g) 
groups <- membership(comm) 
cat("\nDetected groups and their sizes:\n") 
print(table(groups)) 
# Show which documents belong to each group 
groups <- membership(comm) 
split(names(groups), groups) 
plot(comm, g, main = "Communities Detected Using the Louvain Method")

# ============================================================================== 
# ENHANCED DOCUMENT SIMILARITY NETWORK VISUALIZATION 
# ============================================================================== 
                                         
# Attach sentiment scores to network vertices 
V(g)$sentiment <- SentimentA$SentimentQDAP 
# Set vertex size based on eigenvector centrality 
eig_range <- max(eig) - min(eig) 
if (eig_range == 0) { 
  V(g)$size <- rep(15, vcount(g))  # Fixed size if no variance 
} else { 
  V(g)$size <- 12 + 20 * (eig - min(eig)) / eig_range 
} 
# Set vertex color based on sentiment (blue = positive, red = negative) 
sentiment_scaled <- scale(V(g)$sentiment, center = TRUE, scale = TRUE) 
n_colors <- 100 
sentiment_palette <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(n_colors) 
sentiment_bins <- cut(as.numeric(sentiment_scaled), breaks = n_colors,  
                      labels = FALSE, include.lowest = TRUE) 
V(g)$color <- sentiment_palette[sentiment_bins] 
# Set edge width based on connection strength 
E(g)$width <- 1 + E(g)$weight / max(E(g)$weight) * 4 
# Plot enhanced network 
set.seed(33270961) 
plot( 
  g, 
  layout = layout_with_fr, 
  vertex.label.cex = 0.7, 
  vertex.label.color = "darkgreen", 
  vertex.color = V(g)$color, 
  vertex.size = V(g)$size, 
  edge.width = E(g)$width, 
  main = "Improved Single-Mode Network\nColor = SentimentQDAP, Size = Eigenvector, 
Edge Width = Shared Terms" 
) 
# Add sentiment color legend 
legend("bottomleft", 
       legend = c("Negative", "Neutral", "Positive"), 
       fill = sentiment_palette[c(1, n_colors %/% 2, n_colors)], 
       title = "Sentiment") 
                                         
# ============================================================================== 
# TOKEN CO-OCCURRENCE NETWORK (SINGLE-MODE) 
# ============================================================================== 
                                         
# Convert to binary for token analysis 
dtm_bin_token <- (dtm_mat > 0) + 0 
# Compute token-by-token co-occurrence matrix 
adj_mat_token <- t(dtm_bin_token) %*% dtm_bin_token 
# Remove self-loops 
diag(adj_mat_token) <- 0 
# Create igraph object for token network 
g <- graph_from_adjacency_matrix(adj_mat_token, mode = "undirected", weighted = TRUE) 
# Plot basic token co-occurrence network 
set.seed(33270961) 
plot(g, vertex.label.cex = 0.7, main = "Single-Mode Network: Token Co-Occurrence") 
# Calculate centrality measures for tokens 
inv_weights <- 1 / E(g)$weight 
# Degree centrality for tokens 
deg <- degree(g) 
cat("Top degree centrality tokens:\n") 
print(sort(deg, decreasing = TRUE)) 
# Betweenness centrality for tokens 
betw <- betweenness(g, weights = inv_weights) 
cat("\nTop betweenness centrality tokens:\n") 
print(sort(betw, decreasing = TRUE)) 
# Closeness centrality for tokens 
close <- closeness(g, weights = inv_weights) 
cat("\nTop closeness centrality tokens:\n") 
print(sort(close, decreasing = TRUE)) 
# Eigenvector centrality for tokens 
eig <- eigen_centrality(g)$vector 
cat("\nTop eigenvector centrality tokens:\n") 
print(sort(eig, decreasing = TRUE)) 
# Community detection for tokens 
comm <- cluster_louvain(g) 
groups <- membership(comm) 
cat("\nDetected token groups and their sizes:\n") 
print(table(groups)) 
# Show which tokens belong to each group 
groups <- membership(comm) 
split(names(groups), groups) 
plot(comm, g, main = "Token Communities Detected Using the Louvain Method")

# ============================================================================== 
# ENHANCED TOKEN CO-OCCURRENCE NETWORK VISUALIZATION 
# ============================================================================== 
                                         
# Set vertex size based on betweenness centrality 
betw_range_token <- max(betw) - min(betw) 
if (betw_range_token == 0) { 
  V(g)$size <- rep(15, vcount(g)) 
} else { 
  V(g)$size <- 12 + 20 * (betw - min(betw)) / betw_range_token 
} 
# Set vertex color based on closeness centrality (white to dark red) 
n_colors <- 100 
token_palette <- colorRampPalette(c("white", "red4"))(n_colors) 
closeness_bins <- cut(as.numeric(close), breaks = n_colors,  
                      labels = FALSE, include.lowest = TRUE) 
V(g)$color <- token_palette[closeness_bins] 
# Set edge width based on co-occurrence frequency 
E(g)$width <- 1 + E(g)$weight / max(E(g)$weight) * 4 
# Plot enhanced token network 
set.seed(33270961) 
plot( 
  g, 
  layout = layout_with_fr, 
  vertex.label.cex = 0.7, 
  vertex.label.color = "black", 
  vertex.color = V(g)$color, 
  vertex.size = V(g)$size, 
  edge.width = E(g)$width, 
  main = "Improved Single-Mode Network\n(Size = Betweenness, Color = Closeness, Edge 
Width = Token Co-Occurrence Frequency)" 
) 
# Add closeness centrality color legend 
legend("topright", 
       legend = c("Low Closeness", "Medium Closeness", "High Closeness"), 
       fill = token_palette[c(1, n_colors %/% 2, n_colors)], 
       title = "Closeness Centrality")

# ============================================================================== 
# BIPARTITE NETWORK CONSTRUCTION 
# ============================================================================== 
                                         
# Prepare edge list from document-term matrix 
dtmsa <- as.data.frame(dtm_mat) 
dtmsa$doc_id <- rownames(dtmsa) 
dtmsb <- data.frame() 
# Create edge list: iterate through documents and tokens 
for (i in 1:nrow(dtmsa)) { 
  for (j in 1:(ncol(dtmsa) - 1)) {  # Skip doc_id column 
    touse <- cbind( 
      dtmsa[i, j],                  # Weight (term frequency/presence) 
      dtmsa[i, "doc_id"],           # Document ID 
      colnames(dtmsa)[j]            # Token 
    ) 
    dtmsb <- rbind(dtmsb, touse) 
  } 
} 
# Set column names and convert weights to numeric 
colnames(dtmsb) <- c("weight", "doc_id", "token") 
dtmsb$weight <- as.numeric(as.character(dtmsb$weight)) 
# Keep only edges with nonzero weight and reorder columns 
dtmsc <- dtmsb[dtmsb$weight != 0, ] 
dtmsc <- dtmsc[, c(2, 3, 1)]  # Reorder: doc_id, token, weight 
# Create bipartite graph object 
g <- graph_from_data_frame(dtmsc, directed = FALSE) 
V(g)$type <- bipartite_mapping(g)$type 
# Set initial colors and shapes (documents = pink circles, tokens = green squares) 
V(g)$color <- ifelse(V(g)$type, "lightgreen", "pink") 
V(g)$shape <- ifelse(V(g)$type, "circle", "square") 
# Plot basic bipartite network 
set.seed(33270961) 
plot(g, vertex.label.cex = 0.7, main = "Bipartite Document-Token Network")

# ============================================================================== 
# BIPARTITE NETWORK ANALYSIS 
# ============================================================================== 
                                         
# Separate documents and tokens for analysis 
docs_bp <- V(g)[V(g)$type == FALSE]$name  # Documents (type = FALSE) 
tokens_bp <- V(g)[V(g)$type == TRUE]$name  # Tokens (type = TRUE) 
# Calculate degree centrality 
deg <- degree(g) 
cat("\nTop degree centrality - documents:\n") 
print(sort(deg[docs_bp], decreasing = TRUE)) 
cat("\nTop degree centrality - tokens:\n") 
print(sort(deg[tokens_bp], decreasing = TRUE)) 
# Community detection for bipartite network 
comm <- cluster_louvain(g) 
groups <- membership(comm) 
cat("\nDetected bipartite groups and their sizes:\n") 
print(table(groups)) 
plot(comm, g, main = "Bipartite Communities Detected Using the Louvain Method") 
# List nodes in each community 
group_list <- split(names(groups), groups) 
print(group_list)

# ============================================================================== 
# ENHANCED BIPARTITE NETWORK VISUALIZATION 
# ============================================================================== 
                                         
# Add genre colors to documents based on filename patterns 
genres <- ifelse(grepl("politics", V(g)$name), "Politics", 
                 ifelse(grepl("review", V(g)$name), "Review", 
                        ifelse(grepl("sport", V(g)$name), "Sport", NA))) 
# Set color palette for document genres 
doc_palette <- c("Politics" = "skyblue", "Review" = "orange", "Sport" = "purple") 
V(g)$color <- ifelse(V(g)$type, "lightgray", doc_palette[genres]) 
# Set vertex sizes: tokens vary by degree, documents have fixed size 
deg_range_token <- max(deg[tokens_bp]) - min(deg[tokens_bp]) 
if (deg_range_token == 0) { 
  V(g)$size[V(g)$type == TRUE] <- 15  # Fixed size if no variance 
} else { 
  V(g)$size[V(g)$type == TRUE] <- 12 + 8 * (deg[tokens_bp] - min(deg[tokens_bp])) / 
    deg_range_token 
} 
# Set constant size for document nodes 
V(g)$size[V(g)$type == FALSE] <- 15 
# Set edge width based on connection weight 
E(g)$width <- 1 + (E(g)$weight / max(E(g)$weight)) * 8 
# Create vertex labels (only show token names for clarity) 
vertex_labels <- ifelse(V(g)$type, V(g)$name, NA) 
# Plot final enhanced bipartite network 
set.seed(33270961) 
plot( 
  g, 
  layout = layout_with_fr, 
  vertex.label.cex = 0.7, 
  vertex.label = vertex_labels, 
  vertex.label.color = "blue", 
  vertex.color = V(g)$color, 
  vertex.size = V(g)$size, 
  vertex.shape = V(g)$shape, 
  edge.width = E(g)$width, 
  main = "Improved Bipartite Document-Token Network\nColor = Genre, Shape = Node Type, 
Token Node Size = Degree" 
) 
# Add legends for interpretation 
legend("topright",  
       legend = names(doc_palette),  
       fill = doc_palette,  
       title = "Document Genre") 
legend("bottomleft",  
       legend = c("Token", "Document"),  
       pch = c(22, 21),  
       pt.cex = 1.5, 
       col = c("black", "black"),  
       pt.bg = c("white", "white"),  
       title = "Node Type")

# ============================================================================== 
# END OF SCRIPT 
# ==============================================================================
