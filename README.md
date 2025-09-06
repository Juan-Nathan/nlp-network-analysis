# 📚 Natural Language Processing & Network Analysis of a Document Corpus

This project applies natural language processing (NLP) and network analysis techniques in R to analyze a collection of documents from three distinct genres: **Politics**, **Sports**, and **Reviews**.

The workflow involves preprocessing raw text data, generating a document-term matrix, clustering documents, analyzing sentiment, and constructing multiple types of networks—document, token, and bipartite—to explore deeper relationships in the corpus.

## Corpus Overview

- **Total Documents**: 24  
- **Genres**: Politics (8), Sports (8), Reviews (8)  
- **Minimum Document Length**: 200 words  
- **Format**: Plain text files in the `corpus` folder  
- **Sources**:  
  - Politics: e.g., BBC News  
  - Sports: e.g., ESPN  
  - Reviews: e.g., IGN  
- **Naming Convention**: `genre_id.txt` (e.g., `politics_1.txt`)

## Technologies Used

- **Language**: R
- **IDE**: RStudio
- **Packages**: `tm`, `SnowballC`, `proxy`, `SentimentAnalysis`, `igraph`, `RColorBrewer`

## Methodology

### Preprocessing
Performed using `tm`, `SnowballC`, and regular expressions:
- Removed numbers, punctuation, quotes.
- Converted to lowercase.
- Removed common and stop words.
- Applied stemming.

### Document-Term Matrix (DTM)
- Created using `tm`:
  - Converted the preprocessed corpus into a DTM.
  - Removed sparse terms from the DTM (retaining terms in ≥ 30% of documents).
- Manually selected 30 informative tokens to remain in the matrix.

## Analysis & Results

### Hierarchical Clustering
- **Distance**: Cosine  
- **Linkage**: Ward.D  
- **Clustering Accuracy**: 23/24 (≈96%)  
- Documents grouped almost perfectly by genre, indicating strong separability with minimal overlap.
- `sports_5.txt` was misclassified with Reviews due to overlapping analytical language.

### Sentiment Analysis
- **Tool**: `SentimentAnalysis` 
- **Dictionary**: QDAP  
- **Metrics**:  
  - `SentimentQDAP` (net polarity)  
  - `PositivityQDAP` (positive word proportion)  
- **Findings** (from descriptive statistics and hypothesis testing):
  - `SentimentQDAP`:
    - Sports had the most positive overall polarity.
    - Politics showed the lowest and most variable overall polarity.
  - `PositivityQDAP`:
    - Reviews had the highest median positivity and the smallest range, indicating the most consistent use of positive words. 
    - Politics and Sports had similar median positivity and greater variability than Reviews.
    - Reviews were significantly more positive than Politics.
    - No significant difference in positivity was found between Reviews and Sports, or Sports and Politics.

### Single-Mode Document Network
- Nodes = Documents  
- Edges = Number of shared tokens between documents   
- **Important Documents**: `politics_5.txt`, `reviews_3.txt`, `reviews_7.txt`  
- **Communities**:  
  - Documents mainly grouped by genre.
  - Exceptions like `reviews_1.txt`, which grouped with Politics, reflected shared themes or vocabulary.
- **Enhanced Network**: 
  - Node color = `SentimentQDAP`
  - Node size = Eigenvector centrality 
  - Edge width = Shared token count

### Token Co-Occurrence Network
- Nodes = Tokens  
- Edges = Co-occurrence frequency across documents  
- **Important Tokens**: `world`, `fight`, `stori`, `state`
- **Communities**: 
  - Tokens largely grouped by genre.
  - Exceptions like `futur` and `kill`, which grouped with Reviews instead of Politics, reflected overlap in usage across different genres.
- **Enhanced Network**: 
  - Node color = Closeness centrality
  - Node size = Betweenness centrality
  - Edge width = Co-occurrence frequency

### Bipartite Document-Token Network
- Documents linked to tokens they contain
- Nodes = Documents and tokens
- Edges = Token frequency in document
- **Findings**:  
  - Documents and tokens generally grouped by genre.
  - Exceptions like `sports_5.txt` and `event`, which grouped with Reviews, reflected shared themes or vocabulary.
- **Enhanced Network**: 
  - Node color = Genre
  - Node shape = Node type
  - Token node size = Degree
  - Edge width = Token frequency

## Summary

### Important Documents and Tokens
Most documents and tokens were highly interconnected due to shared vocabulary. However, centrality analysis highlighted several important nodes:
- **Documents**: `politics_5.txt`, `reviews_3.txt`, and `reviews_7.txt` consistently showed high centrality, acting as bridges between genres.
- **Tokens**: `world`, `fight`, `stori`, and `state` frequently co-occurred and played key connective roles across the corpus.

### Groups and Clusters
Community detection and hierarchical clustering both grouped documents and tokens primarily by genre. Politics, Reviews, and Sports formed distinct clusters, with few overlaps like `sports_5.txt`, which occasionally grouped with Reviews due to thematic similarity.

### Clustering vs. Network Analysis
- **Clustering** was highly accurate (≈96%) and effective for distinguishing major genre divisions.
- **Network analysis** offered deeper insight into overlapping language, node influence, and structural roles that clustering could not reveal.
- Used together, both methods provide a comprehensive understanding of both dominant groupings and nuanced relationships within the corpus.

## Suggested Improvements

To enhance accuracy and insight:
- Use **lemmatization** over stemming.
- Apply **TF-IDF weighting**.
- Include **n-grams** and **named entities**.
- Use **contextual embeddings** (e.g., BERT).
- Incorporate **POS filtering**, **topic modeling**, **dimensionality reduction**.

## How to Run

1. Clone the repository or download the ZIP file from GitHub.
2. Open the project folder in RStudio.
3. Run the R script (`nlp_network_analysis.r`) inside the RStudio environment.

## Author

Developed by Juan Nathan.
















