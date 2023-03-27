function hello(x, y)
    function no(t)
        return x + y + t
    end
return no
end 

hi = hello(1,2)
display(hi(2))