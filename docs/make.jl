using Fimbul
using JutulDarcy
using Jutul
using Literate
using Documenter
using OrderedCollections

using DocumenterCitations
using DocumenterVitepress
##
cd(@__DIR__)

function dir_to_doc_name(x::String)
    x = replace(x, "_" => " ")
    x = uppercase(x[1:1])*x[2:end]
    return x
end

function get_example_paths(; check_empty = true)
    basepth = joinpath(@__DIR__, "..", "examples")
    examples = OrderedDict()
    examples["storage"] = []
    examples["production"] = []
    examples["validation"] = []
    for excat in readdir(basepth)
        if isdir(joinpath(basepth, excat))
            for exfile in readdir(joinpath(basepth, excat))
                if endswith(exfile, ".jl")
                    if !haskey(examples, excat)
                        examples[excat] = []
                    end
                    filename = first(splitext(exfile))
                    push!(examples[excat], filename)
                end
            end
        end
    end
    if check_empty
        for (k, v) in pairs(examples)
            @assert length(examples[k]) > 0 "No examples found for category $k"
        end
    end
    return examples
end

function parse_tags(text)
    m = match(r"<tags:\s*([^>]+)>", text)
    if !isnothing(m)
        return strip.(split(m.captures[1], ","))
    end
    return nothing
end

function example_path_jl(cname, pth)
    fimbul_dir = realpath(joinpath(@__DIR__, ".."))
    return joinpath(fimbul_dir, "examples", cname, "$pth.jl")
end

function example_path_md(cname, pth)
    return joinpath(@__DIR__, "src", "examples", cname, "$pth.md")
end

function tags_from_file(pth)
    lines = readlines(pth)
    for line in lines
        t = parse_tags(line)
        if !isnothing(t)
            return t
        end
    end
    return nothing
end

function normalize_tag_name(tag::String)
    t = lowercase(strip(tag))
    if t == "production"
        # Keep compatibility with the requested spelling.
        return "proudction"
    end
    return t
end

function all_tags()
    descr = OrderedDict{String, String}()
    descr["validation"] = "Validation examples that compare Fimbul simulations with analytical or benchmark reference cases."
    descr["storage"] = "Examples focused on thermal energy storage systems and operating strategies."
    descr["proudction"] = "Examples focused on geothermal heat production workflows and system performance."
    descr["ates"] = "Aquifer Thermal Energy Storage (ATES) examples."
    descr["btes"] = "Borehole Thermal Energy Storage (BTES) examples."
    descr["ftes"] = "Fracture/Fault Thermal Energy Storage (FTES) examples."
    descr["egs"] = "Enhanced Geothermal System (EGS) examples."
    descr["ags"] = "Advanced Geothermal System (AGS) examples."

    palette = [
        "rgb(228, 26, 28)",
        "rgb(55, 126, 184)",
        "rgb(77, 175, 74)",
        "rgb(152, 78, 163)",
        "rgb(255, 127, 0)",
        "rgb(255, 255, 51)",
        "rgb(166, 86, 40)",
        "rgb(247, 129, 191)"
    ]

    out = OrderedDict{String, NamedTuple{(:desc, :color), Tuple{String, String}}}()
    i = 1
    for (k, v) in pairs(descr)
        out[k] = (desc = v, color = palette[mod1(i, length(palette))])
        i += 1
    end
    return out
end

function infer_tags(category::AbstractString, exname::AbstractString)
    tags = String[]
    cat = lowercase(category)
    name = lowercase(exname)

    if cat == "validation"
        push!(tags, "validation")
    elseif cat == "storage"
        push!(tags, "storage")
    elseif cat == "production"
        push!(tags, "proudction")
    end

    for t in ("ates", "btes", "ftes", "egs", "ags")
        if occursin(t, name)
            push!(tags, t)
        end
    end

    return unique(tags)
end

function tags_for_example(category::AbstractString, exname::AbstractString)
    pth = example_path_jl(category, exname)
    file_tags = tags_from_file(pth)
    tags = isnothing(file_tags) ? infer_tags(category, exname) : map(normalize_tag_name, file_tags)
    return unique(tags)
end

function tag_str(tag::AbstractString)
    return tag_str([tag])
end

function tag_str(tag_names::AbstractVector)
    tags = all_tags()
    s = "``` @raw html\n"
    for tag in tag_names
        t = normalize_tag_name(tag)
        @assert haskey(tags, t) "Unknown tag: $t"
        info = tags[t]
        s *= "<ExampleTag text=\"$t\" color=\"$(info.color)\" />\n"
    end
    s *= "```\n"
    return s
end

function collect_examples_by_tag(; check_empty = false)
    ex_paths = get_example_paths(check_empty = check_empty)
    out = OrderedDict{String, Vector{Tuple{String, String}}}()
    for key in keys(all_tags())
        out[key] = Tuple{String, String}[]
    end
    for (category, example_set) in pairs(ex_paths)
        for exname in example_set
            isfile(example_path_md(category, exname)) || continue
            extags = tags_for_example(category, exname)
            for tag in extags
                @assert haskey(out, tag) "Example $exname in $category has unknown tag $tag"
                push!(out[tag], (exname, category))
            end
        end
    end
    return out
end

function write_tags()
    tags = all_tags()
    outdir = joinpath(@__DIR__, "src", "examples", "overview")
    mkpath(outdir)
    outpth = joinpath(outdir, "example_overview.md")
    ex_tags = collect_examples_by_tag(check_empty = false)
    open(outpth, "w") do io
        println(io, "# Example overview\n")
        println(io, "Fimbul.jl examples are categorized by tags. Use the overview below to quickly find relevant examples by topic and system type.\n")
        for (tag, info) in pairs(tags)
            println(io, "## $tag\n")
            println(io, tag_str(tag))
            println(io, "$(info.desc)\n")
            println(io, "### Examples with the $(lowercase(tag)) tag:\n")
            if length(ex_tags[tag]) == 0
                println(io, "_No examples with this tag yet._\n")
            else
                for (exname, category) in ex_tags[tag]
                    exlink = joinpath("..", "..", "examples", category, "$exname.md")
                    println(io, "1. [$exname]($exlink) (in $category)")
                end
                println(io, "\n")
            end
        end
    end
    println("Wrote tags to $outpth")
end

function timer_str()
    start = "example_t_start = time_ns(); # hide\n"
    stop_1 = "\nt_s = (time_ns() - example_t_start) / 1e9 # hide\n"
    stop = stop_1*"println(\"This example took "*raw"$t_s"*" seconds to complete.\") # hide"
    return (start, stop)
end

function post_run_variables_gc()
    # Attempt of post-processing to GC some objects in temporary @example modules
    # https://discourse.julialang.org/t/delete-a-module/62226
    s =  "\nfunction clear_module!(M::Module)        # hide\n"*
    "    for name ∈ names(M, all=true)        # hide\n"*
    "        if !isconst(M, name)             # hide\n"*
    raw"            @eval M $name = $nothing     # hide"*
    "\n        end                              # hide\n"*
    "    end                                  # hide\n"*
    "end                                      # hide\n"*
    "clear_module!(@__MODULE__)               # hide\n"*
    "GC.gc();                                 # hide\n"
    return s
end

function example_info_footer(subdir, exname)
    return "\n\n# ## Example on GitHub\n"*
    "# If you would like to run this example yourself, it can be downloaded from "*
    "the Fimbul.jl GitHub repository [as a script](https://github.com/sintefmath/Fimbul.jl/blob/main/examples/$subdir/$exname.jl)."
end

function update_footer(content, subdir, exname)
    info_footer = example_info_footer(subdir, exname)
    gc_footer = post_run_variables_gc()
    start, stop = timer_str()
    new_content = string(start, content, info_footer, stop, gc_footer)
    # print(new_content)
    return new_content
end

function replace_tags(content, subdir, exname)
    tags = tags_for_example(subdir, exname)
    if isempty(tags)
        return content
    end
    return string(tag_str(tags), "\n", content)
end

function build_fimbul_docs(
        build_format = nothing;
        build_examples = true,
        build_docs = true,
        build_validation_examples = build_examples,
        build_notebooks = true,
        examples_explicit_list = missing,
        skip_examples = String[],
        clean = true,
        deploy = true,
        use_vitepress = !Sys.iswindows()
    )
    if examples_explicit_list isa String
        examples_explicit_list = [examples_explicit_list]
    end
    examples_explicit_list::Union{Vector{String}, Missing}
    has_explicit_list = !ismissing(examples_explicit_list)
    if has_explicit_list
        @info "Building only examples as examples_explicit_list was specified" examples_explicit_list
    end
    DocMeta.setdocmeta!(Fimbul, :DocTestSetup, :(using Fimbul); recursive=true)
    # DocMeta.setdocmeta!(JutulDarcy, :DocTestSetup, :(using JutulDarcy); recursive=true)
    # DocMeta.setdocmeta!(Jutul, :DocTestSetup, :(using Jutul); recursive=true)
    bib = CitationBibliography(joinpath(@__DIR__, "src", "refs.bib"))

    ## Literate pass
    # Base directory
    fimbul_dir = realpath(joinpath(@__DIR__, ".."))
    # Convert examples as .jl files to markdown
    examples = get_example_paths(check_empty = !has_explicit_list)
    validation_markdown = []
    examples_by_name = OrderedDict{String, Any}()
    if clean
        for (category, example_set) in pairs(examples)
            for ex in example_set
                delpath = joinpath(@__DIR__, "src", "examples", category, "$ex.md")
                if isfile(delpath)
                    println("Deleting generated example \"$ex\":\n\t$delpath")
                    rm(delpath)
                else
                    println("Did not find generated example \"$ex\", skipping removal:\n\t$delpath")
                end
            end
        end
    end
    example_path(cname, pth) = joinpath(fimbul_dir, "examples", cname, "$pth.jl")
    out_dir = joinpath(@__DIR__, "src", "examples")
    notebook_dir = joinpath(@__DIR__, "assets")
    for (category, example_set) in pairs(examples)
        if category == "validation"
            ex_dest = validation_markdown
            do_build = build_validation_examples
        else
            ex_dest = []
            examples_by_name[category] = ex_dest
            do_build = build_examples
        end
        for exname in example_set
            if has_explicit_list
                if exname in examples_explicit_list
                    jutul_message("Examples", "$category/$exname added to build from explicit list.", color = :green)
                else
                    jutul_message("Examples", "$category/$exname not in explicit list, skipping.", color = :yellow)
                    continue
                end
            elseif do_build && !(exname in skip_examples)
                jutul_message("Examples", "$category/$exname was added.", color = :green)
            else
                jutul_message("Examples", "$category/$exname was skipped.", color = :blue)
                continue
            end
            in_pth = example_path(category, exname)
            push!(ex_dest, joinpath("examples", category, "$exname.md"))
            upd(content) = update_footer(content, category, exname)
            fixt(content) = replace_tags(content, category, exname)
            Literate.markdown(in_pth, joinpath(out_dir, category), preprocess = upd, postprocess = fixt)
        end
    end
    examples_markdown = Any["examples/overview/example_overview.md"]
    for (k, v) in pairs(examples_by_name)
        push!(examples_markdown, dir_to_doc_name(k) => v)
    end

    # Must run before makedocs so Documenter validates links against the
    # current set of generated example pages.
    write_tags()

    ## Docs
    if isnothing(build_format)
        println("Use vitepress? ", use_vitepress)
        if use_vitepress
            build_format = DocumenterVitepress.MarkdownVitepress(
                repo = "https://github.com/sintefmath/Fimbul.jl",
            )
        else
            build_format = Documenter.HTML(;
                prettyurls=get(ENV, "CI", "false") == "true",
                canonical="https://sintefmath.github.io/Fimbul.jl",
                edit_link="main",
                size_threshold_ignore = [
                    "ref/jutul.md",
                    "docstrings.md",
                    "man/first_ex.md"
                ],
                assets=String["assets/citations.css"],
            )
        end
    end
    build_pages = [
        "Manual" => [
            "Introduction" => [
                "Fimbul.jl" => "index.md",
            ],
            "Formulation" => [
                "Fluid properties" => "man/formulation/fluid_properties.md",
            ],
            "Cases" => [
                "man/cases/cases.md",
                "man/cases/utils.md",
            ],
            "References" => [
                "Bibliography" => "extras/refs.md"
            ],
        ],
        "Examples" => examples_markdown,
        "Validation" => [
            "Models" => validation_markdown,
        ]
    ]
    # for (k, subpages) in build_pages
    #     println("$k")
    #     @info "$k" subpages
    # end
    if build_docs
        makedocs(;
            modules = [Fimbul],
            authors = "Øystein Klemetsdal <oystein.klemetsdal@sintef.no> and contributors",
            repo = "https://github.com/sintefmath/Fimbul.jl/blob/{commit}{path}#{line}",
            warnonly = false,
            sitename = "Fimbul.jl",
            checkdocs = :exports,
            plugins = [bib],
            format = build_format,
            pages = build_pages,
            draft = get(ENV, "FIMBUL_DOCS_DRAFT_MODE", "0") == "1"
        )
    end
    if build_notebooks
        # Subfolder of final site build folder
        notebook_dir = joinpath(@__DIR__, "build", "final_site", "notebooks")
        mkpath(notebook_dir)
        for (category, example_set) in pairs(examples)
            for exname in example_set
                in_pth = example_path(category, exname)
                ex_notebook_dir = joinpath(notebook_dir, category)
                @info "$exname Writing notebook to $ex_notebook_dir"
                Literate.notebook(in_pth, ex_notebook_dir, execute = false)
            end
        end
    end

    if deploy
        DocumenterVitepress.deploydocs(;
            repo="github.com/sintefmath/Fimbul.jl.git",
            devbranch="main",
            target = "build", # this is where Vitepress stores its output
            branch = "gh-pages",
            push_preview = true
        )
    end
end
# To preview, go to the docs folder and run:
# # DocumenterVitepress.dev_docs("build")
# To only build some examples you can set
# ENV["FIMBUL_DOCS_EXAMPLES_SKIP"] = 1
# You can also enable build after (Linux only):
# ENV["FIMBUL_RUN_VITEPRESS"] = 1
# To use the standard documenter draft mode you can set
# ENV["FIMBUL_DOCS_DRAFT_MODE"] = 1
# This skips all examples, including the inline ones.
if get(ENV, "FIMBUL_DOCS_EXAMPLES_SKIP", "0") == "1"
    # You can add a list of examples to build by running
    # examples_to_build = ["geothermal_1well"]
    if isdefined(Main, :examples_to_build)
        ex_list = examples_to_build
    else
        ex_list = missing
    end
    build_fimbul_docs(
        build_examples = false,
        build_validation_examples = false,
        build_notebooks = false,
        examples_explicit_list = ex_list,
        deploy = false
    )
else
    build_fimbul_docs()
end

if get(ENV, "FIMBUL_RUN_VITEPRESS", "0") == "1" && !Sys.iswindows()
    DocumenterVitepress.dev_docs("build")
end
