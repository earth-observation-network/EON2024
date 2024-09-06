# EON Summer School 2024

Please do not push new/updated `.Rmd` or .`qmd` files that run any code directly, but rather render the website locally (on your own computer), and then push the new files to GitHub.
To render the website locally, you need to use a terminal (not an R console):

```
quarto render
```

This will add new files (pre-rendered versions of the `.Rmd`/`.qmd` documents) to the `_freeze` folder.
You can then push these files to GitHub.
This push will automatically trigger a new GitHub Actions build, which will update the website.
