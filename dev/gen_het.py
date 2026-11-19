# Generates the .a/.r/.u.yaml files of the analyses of hetfamily. Run from the module root:  python dev/gen_het.py
NL = chr(10)


def q(s):
    return "'" + s.replace("'", "''") + "'"


def var(name, title, multi=False, perm=("numeric",), sug=None):
    s = f"    - name: {name}{NL}      title: {q(title)}{NL}      type: {'Variables' if multi else 'Variable'}{NL}"
    if sug:
        s += f"      suggested:{NL}" + "".join(f"          - {x}{NL}" for x in sug)
    s += f"      permitted:{NL}" + "".join(f"          - {x}{NL}" for x in perm) + NL
    return s


def num(name, title, default, mn=None, mx=None, typ="Number"):
    s = f"    - name: {name}{NL}      title: {q(title)}{NL}      type: {typ}{NL}      default: {default}{NL}"
    if mn is not None:
        s += f"      min: {mn}{NL}"
    if mx is not None:
        s += f"      max: {mx}{NL}"
    return s + NL


def txt(name, title, default=""):
    return f"    - name: {name}{NL}      title: {q(title)}{NL}      type: String{NL}      default: {q(default)}{NL}{NL}"


def lst(name, title, opts, default):
    s = f"    - name: {name}{NL}      title: {q(title)}{NL}      type: List{NL}      options:{NL}"
    for n, t in opts:
        s += f"          - name: {n}{NL}            title: {q(t)}{NL}"
    return s + f"      default: {default}{NL}{NL}"


def boo(name, title, default="false"):
    return f"    - name: {name}{NL}      title: {q(title)}{NL}      type: Bool{NL}      default: {default}{NL}{NL}"


def head_a(name, title, mtitle, msub):
    return (f"---{NL}name: {name}{NL}title: {q(title)}{NL}jas: '1.2'{NL}version: '1.0.0'{NL}menuGroup: Heterogeneity{NL}"
            f"menuTitle: {q(mtitle)}{NL}menuSubtitle: {q(msub)}{NL}{NL}options:{NL}    - name: data{NL}      type: Data{NL}{NL}")


def tail_a():
    return boo("showHelp", "Method notes") + boo("showGuide", "Data and setup guide") + "..." + NL


ARMV = lambda: (var("treat", "Treatment (arms)", False, ("numeric", "factor", "id"), ("nominal",)) + txt("control", "Control group value (blank = the first)") + txt("armTreat", "Treatment group value (blank = the second)"))
COVS = lambda: var("covs", "Covariates", True, ("numeric", "factor"), ("continuous", "nominal"))
PROP = lambda: var("prop", "Known propensity score (optional)", False, ("numeric",), ("continuous",))
LEARN = lambda: "".join(boo(n, t, d) for n, t, d in (("lrnLasso", "Lasso", "true"), ("lrnRidge", "Ridge", "false"), ("lrnEnet", "Elastic net", "false"), ("lrnRf", "Random forest", "false"), ("lrnOls", "Linear regression", "false")))
LEARN_NAMES = "lrnLasso, lrnRidge, lrnEnet, lrnRf, lrnOls"
RULE = lambda: lst("lambdaRule", "Lasso and ridge penalty", [("min", "Cross-validation minimum"), ("1se", "One standard error above the minimum")], "min")


def gen_a(key):
    if key == "hetgeneric":
        s = head_a("hetgeneric", "Heterogeneous Effects (BLP, GATES, CLAN)", "Heterogeneous Effects", "Best linear predictor, sorted group effects and classification analysis")
        s += var("outcome", "Outcome") + ARMV() + COVS() + PROP() + LEARN()
        s += num("groups", "Number of groups (GATES, CLAN)", 5, 2, 10, "Integer") + num("nsplit", "Number of sample splits", 50, 5, 500, "Integer") + num("shareMain", "Share of the sample used for estimation", 0.5, 0.2, 0.8)
        s += num("level", "Confidence level (%)", 95, 50, 99.9) + lst("vcov", "Standard errors", [("hc3", "Robust (HC3)"), ("hc1", "Robust (HC1)"), ("const", "Classical")], "hc3") + RULE()
        s += boo("monotonize", "Sort the group effects (monotonize)", "true") + boo("clan", "Classification analysis", "true") + boo("plot", "Group effects plot", "true")
        s += num("seed", "Random seed", 1234, 1, 999999, "Integer")
    elif key == "hetprofile":
        s = head_a("hetprofile", "Effect Profile (CATE regression)", "Effect Profile", "How the effect changes with each covariate")
        s += var("outcome", "Outcome") + ARMV() + COVS() + var("moderators", "Moderators (default: all covariates)", True, ("numeric", "factor"), ("continuous", "nominal")) + PROP() + LEARN()
        s += num("nsplit", "Number of sample splits", 50, 5, 500, "Integer") + num("shareMain", "Share of the sample used for estimation", 0.5, 0.2, 0.8)
        s += num("level", "Confidence level (%)", 95, 50, 99.9) + lst("vcov", "Standard errors", [("hc3", "Robust (HC3)"), ("hc1", "Robust (HC1)"), ("const", "Classical")], "hc3") + RULE()
        s += boo("plot", "Coefficient plot", "true") + num("seed", "Random seed", 1234, 1, 999999, "Integer")
    elif key == "hetrate":
        s = head_a("hetrate", "Priority Curve (RATE)", "Priority Curve", "Does ranking by predicted effect put the largest gains first?")
        s += var("outcome", "Outcome") + ARMV() + COVS() + PROP() + LEARN()
        s += num("nsplit", "Number of sample splits", 50, 5, 500, "Integer") + num("nboot", "Bootstrap samples per split", 100, 20, 2000, "Integer") + num("shareMain", "Share of the sample used for estimation", 0.5, 0.2, 0.8)
        s += num("level", "Confidence level (%)", 95, 50, 99.9) + RULE() + boo("plot", "Priority curve plot", "true") + num("seed", "Random seed", 1234, 1, 999999, "Integer")
    else:
        s = head_a("hetpolicy", "Treatment Rule Value", "Treatment Rule", "Value of treating the units with the largest predicted effect")
        s += var("outcome", "Outcome") + ARMV() + COVS() + PROP() + LEARN()
        s += num("cost", "Cost of treating a unit (outcome units)", 0) + num("nsplit", "Number of sample splits", 50, 5, 500, "Integer") + num("shareMain", "Share of the sample used for estimation", 0.5, 0.2, 0.8)
        s += num("level", "Confidence level (%)", 95, 50, 99.9) + RULE() + num("seed", "Random seed", 1234, 1, 999999, "Integer")
    return s + tail_a()


def col(name, title, typ, extra=""):
    return f"        - name: {name}{NL}          title: {q(title)}{NL}          type: {typ}{NL}" + (extra + NL if extra else "")


def tbl(name, title, clear, cols, visible=None):
    s = f"    - name: {name}{NL}      title: {q(title)}{NL}      type: Table{NL}      rows: 0{NL}"
    if visible:
        s += f"      visible: {visible}{NL}"
    s += f"      clearWith: [{clear}]{NL}      columns:{NL}" + "".join(cols) + NL
    return s


PC = "          format: zto,pvalue"


def gen_r(key, title):
    s = (f"---{NL}name: {key}{NL}title: {q(title)}{NL}jrs: '1.1'{NL}{NL}items:{NL}"
         f"    - name: guide{NL}      title: Data and setup guide{NL}      type: Html{NL}      visible: (showGuide){NL}{NL}"
         f"    - name: helpText{NL}      title: Method notes{NL}      type: Html{NL}      visible: (showHelp){NL}{NL}"
         f"    - name: notes{NL}      title: Notes{NL}      type: Html{NL}      visible: true{NL}{NL}")
    stats = [col("est", "Estimate", "number"), col("se", "Std. error", "number"), col("lo", "CI lower", "number"), col("hi", "CI upper", "number"), col("p", "p", "number", PC)]
    if key == "hetgeneric":
        c = f"outcome, treat, control, armTreat, covs, prop, {LEARN_NAMES}, groups, nsplit, shareMain, level, vcov, lambdaRule, monotonize, seed"
        s += tbl("fit", "How well each learner predicts the effects", c, [col("learner", "Learner", "text"), col("lam", "Lambda (BLP)", "number"), col("lamb", "Lambda-bar (GATES)", "number"),
                                                                        col("bb", "Best for BLP", "text"), col("bg", "Best for GATES and CLAN", "text")])
        s += tbl("blp", "Best linear predictor of the effects", c, [col("what", "", "text")] + stats)
        s += tbl("gates", "Average effect by group of predicted effect (GATES)", c, [col("what", "Group", "text")] + stats)
        s += tbl("clan", "Who is affected most: covariate means (CLAN)", c + ", clan", [col("what", "Covariate", "text"), col("m1", "Least affected", "number"), col("m2", "Most affected", "number"), col("est", "Difference", "number"),
                                                                                        col("se", "Std. error", "number"), col("lo", "CI lower", "number"), col("hi", "CI upper", "number"), col("p", "p", "number", PC)], "(clan)")
        s += (f"    - name: gatesPlot{NL}      title: 'Average effect by group'{NL}      type: Image{NL}      width: 560{NL}      height: 360{NL}"
              f"      renderFun: .plotGates{NL}      visible: (plot){NL}      clearWith: [{c}]{NL}{NL}")
        refs = "cddf2018, genericml"
    elif key == "hetprofile":
        c = f"outcome, treat, control, armTreat, covs, moderators, prop, {LEARN_NAMES}, nsplit, shareMain, level, vcov, lambdaRule, seed"
        s += tbl("fit", "How well each learner predicts the effects", c, [col("learner", "Learner", "text"), col("lam", "Lambda (BLP)", "number"), col("lamb", "Lambda-bar (GATES)", "number"), col("bg", "Used for the profile", "text")])
        s += tbl("prof", "Effect of the treatment and how it changes with each moderator (per standard deviation)", c, [col("what", "", "text")] + stats)
        s += (f"    - name: profPlot{NL}      title: 'Change of the effect per standard deviation of each moderator'{NL}      type: Image{NL}      width: 560{NL}      height: 420{NL}"
              f"      renderFun: .plotProf{NL}      visible: (plot){NL}      clearWith: [{c}]{NL}{NL}")
        refs = "semenova2021, cddf2018"
    elif key == "hetrate":
        c = f"outcome, treat, control, armTreat, covs, prop, {LEARN_NAMES}, nsplit, nboot, shareMain, level, lambdaRule, seed"
        s += tbl("fit", "How well each learner predicts the effects", c, [col("learner", "Learner", "text"), col("lam", "Lambda (BLP)", "number"), col("lamb", "Lambda-bar (GATES)", "number"), col("bg", "Used for the ranking", "text")])
        s += tbl("rate", "Rank-weighted average treatment effects", c, [col("what", "", "text")] + stats)
        s += tbl("toc", "Gain from treating the top share first (TOC)", c, [col("q", "Share treated", "number")] + stats)
        s += (f"    - name: tocPlot{NL}      title: 'Treatment operating characteristic'{NL}      type: Image{NL}      width: 560{NL}      height: 360{NL}"
              f"      renderFun: .plotToc{NL}      visible: (plot){NL}      clearWith: [{c}]{NL}{NL}")
        refs = "yadlowsky2021, cddf2018"
    else:
        c = f"outcome, treat, control, armTreat, covs, prop, {LEARN_NAMES}, cost, nsplit, shareMain, level, lambdaRule, seed"
        s += tbl("fit", "How well each learner predicts the effects", c, [col("learner", "Learner", "text"), col("lam", "Lambda (BLP)", "number"), col("lamb", "Lambda-bar (GATES)", "number"), col("bg", "Used for the rule", "text")])
        s += tbl("val", "Value of the rule: treat a unit if its predicted effect exceeds the cost", c, [col("what", "", "text")] + stats)
        refs = "cddf2018, kitagawatetenov2018"
    return s + f"refs:{NL}" + "".join(f"    - {r.strip()}{NL}" for r in refs.split(",")) + f"{NL}...{NL}"


def box(label, name, single=True):
    s = f"          - type: TargetLayoutBox{NL}            label: {q(label)}{NL}            children:{NL}                - type: VariablesListBox{NL}                  name: {name}{NL}                  isTarget: true{NL}                  height: small{NL}"
    if single:
        s += f"                  maxItemCount: 1{NL}"
    return s


def tb(name, label, fmt="number"):
    return f"          - type: TextBox{NL}            name: {name}{NL}            label: {q(label)}{NL}            format: {fmt}{NL}"


def cb(name, label=None):
    return f"          - type: CheckBox{NL}            name: {name}{NL}" + (f"            label: {q(label)}{NL}" if label else "")


def cmb(name):
    return f"          - type: ComboBox{NL}            name: {name}{NL}"


def collapse(label, body, collapsed=True):
    return f"    - type: CollapseBox{NL}      label: {q(label)}{NL}      collapsed: {str(collapsed).lower()}{NL}      children:{NL}" + body


def gen_u(key, title):
    boxes = [("Outcome", "outcome", 1), ("Treatment (arms)", "treat", 1), ("Covariates", "covs", 0), ("Known propensity score (optional)", "prop", 1)]
    if key == "hetprofile":
        boxes = [("Outcome", "outcome", 1), ("Treatment (arms)", "treat", 1), ("Covariates", "covs", 0), ("Moderators (default: all covariates)", "moderators", 0), ("Known propensity score (optional)", "prop", 1)]
    s = (f"title: {q(title)}{NL}name: {key}{NL}jus: '3.0'{NL}stage: 0{NL}compilerMode: tame{NL}children:{NL}"
         f"    - type: VariableSupplier{NL}      persistentItems: false{NL}      stretchFactor: 1{NL}      children:{NL}")
    for lab, nm, sg in boxes:
        s += box(lab, nm, bool(sg))
    s += collapse("Treatment groups", tb("control", "Control group value", "string") + tb("armTreat", "Treatment group value", "string"))
    s += collapse("Learners", "".join(cb(n) for n in ("lrnLasso", "lrnRidge", "lrnEnet", "lrnRf", "lrnOls")) + cmb("lambdaRule"), False)
    if key == "hetgeneric":
        s += collapse("Estimation", tb("groups", "Number of groups (GATES, CLAN)") + tb("nsplit", "Number of sample splits") + tb("shareMain", "Share of the sample used for estimation") + tb("level", "Confidence level (%)")
                      + cmb("vcov") + cb("monotonize") + tb("seed", "Random seed"))
        s += collapse("Output tables", cb("clan") + cb("plot"))
    elif key == "hetprofile":
        s += collapse("Estimation", tb("nsplit", "Number of sample splits") + tb("shareMain", "Share of the sample used for estimation") + tb("level", "Confidence level (%)") + cmb("vcov") + tb("seed", "Random seed"), False)
        s += collapse("Output tables", cb("plot"))
    elif key == "hetrate":
        s += collapse("Estimation", tb("nsplit", "Number of sample splits") + tb("nboot", "Bootstrap samples per split") + tb("shareMain", "Share of the sample used for estimation") + tb("level", "Confidence level (%)") + tb("seed", "Random seed"), False)
        s += collapse("Output tables", cb("plot"))
    else:
        s += collapse("Estimation", tb("cost", "Cost of treating a unit (outcome units)") + tb("nsplit", "Number of sample splits") + tb("shareMain", "Share of the sample used for estimation") + tb("level", "Confidence level (%)") + tb("seed", "Random seed"), False)
    s += collapse("Output", cb("showHelp", "Method notes") + cb("showGuide", "Data and setup guide"))
    return s


TITLES = {"hetgeneric": "Heterogeneous Effects (BLP, GATES, CLAN)", "hetprofile": "Effect Profile (CATE regression)", "hetrate": "Priority Curve (RATE)", "hetpolicy": "Treatment Rule Value"}
for k, t in TITLES.items():
    for ext, txt_ in (("a", gen_a(k)), ("r", gen_r(k, t)), ("u", gen_u(k, t))):
        with open(f"jamovi/{k}.{ext}.yaml", "w", encoding="utf-8", newline="") as f:
            f.write(txt_)
print("written")
