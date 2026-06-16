# Documentation overview of pseudoknots directory

The directory is managed to have linear dependency.

````mermaid
flowchart TD
    A[parser.rs] --> B[closed_region_tree.rs]
    B --> C[loops.rs]
    C --> D[enumerate.rs]
   

````

