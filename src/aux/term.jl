function xpsmsg(text, title, subtitle; adjust=:left)
    println("\n")
    return Panel(highlight(text);
                 width = 77,
                 justify = adjust,
                 title=title,
                 title_style="red bold",
                 title_justify=:right,
                 subtitle=subtitle,
                 subtitle_style="#9558B2",
                 subtitle_justify=:left)
end

function separator(w = 2)
    # colors are from https://github.com/JuliaLang/julia-logo-graphics
    w = (2 ≤ w ≤ 25) ? w : 2
    gap = 75 - 3w
    println(' '^gap *
        Panel(width=w, box=:HEAVY, style="#389826") *
        Panel(width=w, box=:HEAVY, style="#9558B2") *
        Panel(width=w, box=:HEAVY, style="#CB3C33"))
end
