include("ctRoot.jl")

function replicate_all(;saveall=false)

    fig1 = plot_announcements(slides = false)
    saveall && savefig(fig1, "Graphs/current/announcements_paper.pdf", width = 900, height = 450)

    mt = load("Output/mt.jld2", "mt")

    fig2 = Lplot(mt.ct, slides = false)
    saveall && savefig(fig2, "Graphs/current/first_L_paper.pdf", width = 900, height = 450)

    fig3 = hplot(mt.ct, slides = false)
    saveall && savefig(fig3, "Graphs/current/first_g_paper.pdf", width = 900, height = 350)

    fig4 = Eplot(mt.ct, slides = false)
    saveall && savefig(fig4, "Graphs/current/first_p_paper.pdf", width = 900, height = 400)

    fig5 = plansp(mt, slides = false)
    saveall && savefig(fig5, "Graphs/current/plans_paper.pdf", width = 900, height = 400)

    fig6 = Lωplot(mt, slides = false)
    saveall && savefig(fig6, "Graphs/current/contour_paper.pdf", width = 900, height = 500)

    fig7 = Cplot(mt, slides = false)
    saveall && savefig(fig7, "Graphs/current/Ccontour_paper.pdf", width = 900, height = 450)

    fig8 = avgplans(mt, 100, CIs = false, slides = false)
    saveall && savefig(fig8, "Graphs/current/mimics_paper.pdf", width = 900, height = 400)

    mt = load("Output/mt.jld2", "mt")

    fig9a, fig9b = strategy_μ(mt, slides = false, save_stats = true);
    saveall && savefig(fig9a, "Graphs/current/marg_achi_paper.pdf", width = 600, height = 450)
    saveall && savefig(fig9b, "Graphs/current/marg_omegachi_paper.pdf", width = 600, height = 450)

    # Comparative statics would go here

    fig12 = comp_plot_planner(mt, slides = false)
    saveall && savefig(fig12, "Graphs/current/ramsey_paper.pdf", width = 900, height = 350)

    fig_app1 = Lplot_fixed_ω(mt, slides = false)
    saveall && savefig(fig_app1, "Graphs/current/contour_app_paper.pdf", width = 900, height = 400)


    # Gradual Feedback Rules

    mtp = load("Output/mtpsi_s20_20230720.jld2", "mt_psi");

    fig13a = twolines(mtp, jp = 2, slides = false)
    saveall && savefig(fig13a, "Graphs/current/psithing_L1.pdf", width = 500, height = 350)

    fig13b = twolines(mtp, show = "C", jp = 2, slides = false)
    saveall && savefig(fig13b, "Graphs/current/psithing_C1.pdf", width = 500, height = 350)

    fig14 = implied_plan(mtp, slides = false)
    saveall && savefig(fig14, "Graphs/current/psithing_plans.pdf", width = 1000, height = 450)
end

function thumbnail(;saveall = false)

    mt = load("Output/mt.jld2", "mt")

    fig = Lplot(mt.ct, slides = true, plot_bgcolor="#fafafa", paper_bgcolor="#fafafa")
    saveall && savefig(fig, "cddp_fig.png", width = 850, height = 500)
end
