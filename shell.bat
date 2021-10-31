:loop
    git add . 
    git commit -m"Auto update" 
    git push origin master 
    timeout 100

goto loop