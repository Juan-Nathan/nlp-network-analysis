# 📚 Natural Language Processing & Network Analysis on a Document Corpus

This project applies natural language processing (NLP) and network analysis techniques in R to analyze a collection of documents from three distinct genres: **Politics**, **Sport**, and **Reviews**. The workflow involves preprocessing raw text data, generating a document-term matrix (DTM), clustering documents, analyzing sentiment, and constructing multiple types of networks—document, token, and bipartite—to explore deeper relationships in the corpus.

## Corpus Overview

- **Total Documents**: 24  
- **Genres**: Politics (8), Sport (8), Reviews (8)  
- **Minimum Document Length**: 200 words  
- **Format**: Plain text files in the `corpus` folder  
- **Sources**:  
  - Politics: e.g., BBC News  
  - Sport: e.g., ESPN  
  - Reviews: e.g., IGN  
- **Naming Convention**: `genre_id.txt` (e.g., `politics_1.txt`)

## Technologies Used

- **Language**: R
- **IDE**: RStudio
- **Packages**: `tm`, `SnowballC`, `proxy`, `SentimentAnalysis`, `igraph`, `RColorBrewer`

## Methodology

### Preprocessing
Performed using `tm`, `SnowballC`, and regular expressions:
- Removed numbers, punctuation, quotes
- Converted to lowercase
- Removed custom common words and English stopwords
- Applied stemming

### Document-Term Matrix (DTM)
- Created using `tm`:
  - Converted the preprocessed corpus into a DTM
  - Removed sparse terms from the DTM (retaining terms in ≥ 30% of documents)
- Converted the DTM into a matrix
- Selected 30 informative tokens manually

## Analysis & Results

### Hierarchical Clustering
- **Distance**: Cosine  
- **Linkage**: Ward.D  
- **Clustering Accuracy**: 23/24 (≈96%)  
- Documents grouped almost perfectly by genre, indicating strong separability with minimal overlap
- `sport_5.txt` was misclassified with reviews due to overlapping analytical language

### Sentiment Analysis
- **Tool**: `SentimentAnalysis` 
- **Dictionary**: QDAP  
- **Metrics**:  
  - `SentimentQDAP` (net polarity)  
  - `PositivityQDAP` (positive word proportion)  
- **Findings (from descriptive statistics and hypothesis testing)**:
  - `SentimentQDAP`:
    - Sports articles had the most positive overall tone (highest polarity)
    - Reviews were slightly positive on average
    - Politics had the lowest and most variable polarity scores
  - `PositivityQDAP`:
    - Reviews used positive words most consistently and had the highest median positivity
    - Politics and sports had similar average positivity but showed greater variability
    - Reviews were significantly more positive than politics (p = 0.031)
    - No significant difference between reviews and sports, or sports and politics

### Single-Mode Document Network
- Nodes = Documents  
- Edges = Number of shared tokens between documents   
- **Central Documents**:  
  - `politics_5.txt`, `review_3.txt`, `review_7.txt`  
- **Communities**:  
  - Documents mainly grouped by genre
  - Exceptions like `review_1.txt` reflected shared themes or vocabulary
- **Enhanced Network**: 
  - Node color = sentiment
  - Node size = eigenvector centrality 
  - Edge width = shared token count

### Token Co-Occurrence Network
- Nodes = Tokens  
- Edges = Co-occurrence frequency across documents  
- **Important Tokens**:  
  - `world`, `fight`, `stori`, `state`
- **Communities**: 
  - Tokens largely grouped by genre
  - Exceptions like `futur` and `kill` reflected overlap in usage across different genres
- **Enhanced Network**: 
  - Node color = closeness
  - Node size = betweenness centrality
  - Edge width = co-occurrence frequency

### Bipartite Document-Token Network
- Documents linked to tokens they contain
- Nodes: Documents and tokens
- Edges = Token frequency in document
- **Findings**:  
  - Documents and tokens generally grouped by genre  
  - Exceptions like `sport_5.txt` and `event` reflected shared themes or vocabulary
- **Enhanced Network**: 
  - Node color = genre
  - Node shape = node type
  - Token node size = degree
  - Edge width = token frequency

## Summary

### Important Documents and Tokens
Most documents and tokens were highly interconnected due to shared vocabulary. However, centrality analysis highlighted several key nodes:
- **Documents**: `politics_5.txt`, `review_3.txt`, and `review_7.txt` consistently showed high centrality, acting as bridges between genres.
- **Tokens**: `world`, `fight`, `stori`, and `state` frequently co-occurred and played key connective roles across the corpus.

### Groups and Clusters
Community detection and hierarchical clustering both grouped documents and tokens primarily by genre. Politics, reviews, and sports formed distinct clusters, with few overlaps like `sport_5.txt` appearing in review clusters due to thematic similarity.

### Clustering vs. Network Analysis
- **Clustering** was highly accurate (~96%) and effective for detecting major genre divisions.
- **Network analysis** offered deeper insight into overlapping language, node influence, and structural roles that clustering could not reveal.
- Used together, both methods provide a comprehensive understanding of both dominant groupings and nuanced relationships within the corpus.

## Suggested Improvements

To enhance accuracy and insight:
- Use **lemmatization** over stemming
- Apply **TF-IDF weighting**
- Include **n-grams** and **named entities**
- Use **contextual embeddings** (e.g., BERT)
- Incorporate **POS filtering**, **topic modeling**, **dimensionality reduction**

## How to Run

1. Clone the repository or download the ZIP file from GitHub.
2. Open the project folder in RStudio.
3. Run the R script (`nlp_network_analysis.r`) inside the RStudio environment.

## Author

Developed by Juan Nathan.







