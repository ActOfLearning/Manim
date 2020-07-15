from manimlib.imports import *
import numpy as np

class Cycloid(ParametricFunction):
    CONFIG = {
        "point_a"       : ORIGIN,
        "radius"        : 1,
        "start_theta"   : 0,
        "end_theta"     : 2 * np.pi,
        "density"       : 5 * DEFAULT_POINT_DENSITY_1D,
        "color"         : BLUE,
        "inverted"      : False,
        "frac"          : 1.0
    }
    def __init__(self, **kwargs):
        digest_config(self, kwargs)
        ParametricFunction.__init__(self, self.pos_func, **kwargs)

    def pos_func(self, t):
        T = self.start_theta + t * (self.end_theta - self.start_theta)
        return self.point_a + np.array([
            self.radius * T - (-1 if self.inverted else 1) * self.frac * self.radius * np.sin(T),
            self.radius - self.radius * self.frac * np.cos(T),
            0
        ])

class LearningScene(GraphScene, MovingCameraScene):
    CONFIG = {
        "x_min": -10,
        "x_max": 10,
        "x_axis_width": 20,
        "x_tick_frequency": 1,
        "x_leftmost_tick": None,  # Change if different from x_min
        "x_labeled_nums": None,
        "x_axis_label": "",
        "y_min": -5,
        "y_max": 5,
        "y_axis_height": 10,
        "y_tick_frequency": 1,
        "y_bottom_tick": None,  # Change if different from y_min
        "y_labeled_nums": None,
        "y_axis_label": "",
        "axes_color": BLACK,
        "graph_origin": ORIGIN,
        "exclude_zero_label": True,
        "default_graph_colors": [BLUE, GREEN, YELLOW],
        "default_derivative_color": GREEN,
        "default_input_color": YELLOW,
        "default_riemann_start_color": BLUE,
        "default_riemann_end_color": GREEN,
        "area_opacity": 0.8,
        "num_rects": 50,
        "dot_kwargs": {
            "radius": 0.05,
        },
        "line_kwargs": {
            "stroke_width": 2,
        },
        "fill_triangle_kwargs": {
            # "fill_color": BLUE,
            "fill_opacity": .5,
            "stroke_width": 0,
        },
    }

    def setup(self):
        GraphScene.setup(self)
        MovingCameraScene.setup(self)
    
    def construct(self):
        f = 1
        self.setup_axes()
        self.play(FadeOut(self.axes))
        self.wait()
        self.play(self.camera_frame.shift, UP)
        titl = TexMobject("\\text{\\underline{Cycloid}}").scale(1.5).move_to(4.25 * UP)
        self.play(Write(titl))
        theta = ValueTracker(-2 * PI)
        cyc = Cycloid(frac = f, point_a = -2 * PI * RIGHT, end_theta = 4 * PI)
        b_line = Line(10 * LEFT, 10 * RIGHT).set_color(GREY)
        self.play(Write(b_line))
        self.wait()
        ccirc = Circle()
        cline = Line(ORIGIN, DOWN).rotate_about_origin(-theta.get_value())
        cdot = Dot().move_to(cline.get_end())
        c_grp = VGroup(ccirc, cline, cdot).set_color(RED)
        c_grp.shift(UP + theta.get_value() * RIGHT)
        self.play(Write(c_grp))
        def cupd(obj):
            ncirc = Circle()
            nline = Line(ORIGIN, DOWN).rotate_about_origin(-theta.get_value())
            ndot = Dot().move_to(nline.get_end())
            n_grp = VGroup(ncirc, nline, ndot).set_color(RED)
            n_grp.shift(UP + theta.get_value() * RIGHT)
            obj.become(n_grp)
        c_grp.add_updater(cupd)
        self.add(c_grp)
        self.play(theta.increment_value, 4 * PI, ShowCreation(cyc), run_time = 10)
        self.play(theta.increment_value, -2 * PI, run_time = 10)
        c_grp.clear_updaters()
        self.play(*[FadeOut(i) for i in [b_line, cyc, c_grp]])
        self.play(
            self.camera_frame.move_to, PI * RIGHT + 0.5 * DL + 0.5 * LEFT,
            self.camera_frame.scale, 0.5,
            titl.scale_in_place, 0.5,
            titl.move_to, PI * RIGHT + 0.5 * DL + 1.625 * UP + 0.5 * LEFT,)
        brachis_cyc = Cycloid(end_theta = PI).set_color(BLUE)
        brachis_cyc.rotate_about_origin(PI, axis = RIGHT)
        self.play(Write(brachis_cyc))
        sli_dot = Dot().set_color(RED)
        brachis_title = VGroup(
            TextMobject("Brachistochrone Curve").set_color(YELLOW).scale(0.4375),
            TextMobject("brachis $meaning$ least").scale(0.4375),
            TextMobject("chrone $meaning$ time").scale(0.4375)).arrange(direction = DOWN, buff = 0.125)
        for t in brachis_title[1:]:
            t.align_to(brachis_title[0], LEFT)
        brachis_title.shift(3.5 * RIGHT + 0.5 * DOWN)
        self.play(Write(brachis_title))
        self.play(MoveAlongPath(sli_dot, brachis_cyc), run_time = 5)
        self.play(FadeOut(sli_dot))
        cycl_grp = VGroup(
            Cycloid(end_theta = PI).rotate_about_origin(PI, axis = RIGHT),
            Cycloid(start_theta = PI / 5, end_theta = PI).rotate_about_origin(PI, axis = RIGHT),
            Cycloid(start_theta = 2 * PI / 5, end_theta = PI).rotate_about_origin(PI, axis = RIGHT),
            Cycloid(start_theta = 3 * PI / 5, end_theta = PI).rotate_about_origin(PI, axis = RIGHT),
            Cycloid(start_theta = 4 * PI / 5, end_theta = PI).rotate_about_origin(PI, axis = RIGHT),
        )
        def cycloidpts(t, radius = 1, frac = 1, inverted = False):
            return (radius * t - (-1 if inverted else 1) * frac * radius * np.sin(t)) * RIGHT + (radius - radius * frac * np.cos(t)) * UP
        dot_grp = VGroup()
        for i in range(5):
            temp = Dot().move_to(cycloidpts(i * PI / 5))
            temp.rotate_about_origin(PI, axis = RIGHT)
            dot_grp.add(temp)
        dot_grp[0].set_color(GREEN)
        dot_grp[1].set_color(YELLOW)
        dot_grp[2].set_color(GOLD)
        dot_grp[3].set_color(PURPLE)
        dot_grp[4].set_color(ORANGE)
        self.play(Write(dot_grp))
        tauto_title = VGroup(
            TextMobject("Tautochrone Curve").set_color(YELLOW).scale(0.4375),
            TextMobject("tauto $meaning$ equal").scale(0.4375),
            TextMobject("chrone $meaning$ time").scale(0.4375)).arrange(direction = DOWN, buff = 0.125)
        for t in tauto_title[1:]:
            t.align_to(tauto_title[0], LEFT)
        tauto_title.shift(3.5 * RIGHT + 0.5 * DOWN)
        self.play(Transform(brachis_title, tauto_title))
        self.play(*[MoveAlongPath(d, p) for d, p in zip(dot_grp, cycl_grp)], run_time = 5)
        self.wait()
        self.play(*[FadeOut(i) for i in [brachis_title, titl, dot_grp, brachis_cyc]])
        self.wait(2)

class LearningCava(GraphScene, MovingCameraScene):
    CONFIG = {
        "x_min": -10,
        "x_max": 10,
        "x_axis_width": 20,
        "x_tick_frequency": 1,
        "x_leftmost_tick": None,  # Change if different from x_min
        "x_labeled_nums": None,
        "x_axis_label": "",
        "y_min": -5,
        "y_max": 5,
        "y_axis_height": 10,
        "y_tick_frequency": 1,
        "y_bottom_tick": None,  # Change if different from y_min
        "y_labeled_nums": None,
        "y_axis_label": "",
        "axes_color": GREY,
        "graph_origin": ORIGIN,
        "exclude_zero_label": True,
        "default_graph_colors": [BLUE, GREEN, YELLOW],
        "default_derivative_color": GREEN,
        "default_input_color": YELLOW,
        "default_riemann_start_color": BLUE,
        "default_riemann_end_color": GREEN,
        "area_opacity": 0.8,
        "num_rects": 50,
        "dot_kwargs": {
            "radius": 0.05,
        },
        "line_kwargs": {
            "stroke_width": 2,
        },
        "fill_triangle_kwargs": {
            # "fill_color": BLUE,
            "fill_opacity": .5,
            "stroke_width": 0,
        },
    }

    def setup(self):
        GraphScene.setup(self)
        MovingCameraScene.setup(self)
    
    def construct(self):
        self.setup_axes()
        self.play(FadeOut(self.axes))
        self.wait()
        titl = TexMobject("\\text{\\underline{Cavalieri's Principle}}").move_to(3 * RIGHT + 3.625 * UP)
        self.play(Write(titl), self.camera_frame.move_to, 3 * RIGHT + 1.5 * UP, self.camera_frame.scale, 2 / 3)
        self.wait()
        h = ValueTracker(0)
        base_gr = self.get_graph(lambda t: 0, x_min = 0, x_max = 2)
        hei_gr = self.get_graph(lambda t: 0, x_min = 0, x_max = 2)
        def cava_rectanges(g1, g2, x_min = 0, x_max = 1, dx = 0.1, ashift = 1, stroke_width = 1):
            rectangles = VGroup()
            scolor, ecolor = GREEN, BLUE
            x_range = np.arange(x_min, x_max, dx)
            colors = color_gradient([scolor, ecolor], len(x_range))
            for x, color in zip(x_range, colors):
                lg_pt = self.input_to_graph_point(x, g1)
                rg_pt = self.input_to_graph_point(x, g2) + ashift * UP
                points = VGroup(*list(map(VectorizedPoint, [
                lg_pt, lg_pt + 0.975 * dx * RIGHT, rg_pt])))
                rect = Rectangle()
                rect.replace(points, stretch = True)
                rect.set_fill(color, opacity = 1)
                rect.set_stroke(BLACK, width = stroke_width)
                rectangles.add(rect)
            #rectangles.add(g1.copy())
            #rectangles.add(g2.copy().shift(ashift * UP))
            rectangles.rotate_about_origin(PI, axis = UP)
            rectangles.rotate_about_origin(-PI / 2)
            return rectangles
        def cava_list(g1, g2, n_iterations = 6, max_dx = 0.5, power_base = 2, stroke_width = 1, ashift = 1, x_min = 0, x_max = 1):
            return [cava_rectanges(
                g1 = g1,
                g2 = g2,
                dx = float(max_dx) / power_base ** n,
                stroke_width = float(stroke_width) / power_base ** n,
                ashift = ashift,
                x_min = x_min,
                x_max = x_max
            ) for n in range(n_iterations)]
        def xyreflect(obj):
            res_obj = obj.copy()
            res_obj.rotate_about_origin(PI, axis = UP)
            res_obj.rotate_about_origin(-PI / 2)
            return res_obj
        rect_list = cava_list(base_gr, hei_gr, ashift = 1, x_max = 2, n_iterations = 5)
        rect_list_cpy = rect_list.copy()
        rects = rect_list_cpy[0]
        srect = VMobject()
        srect.set_points_as_corners([ORIGIN, RIGHT, RIGHT + 2 * UP, 2 * UP, ORIGIN]).set_color(BLUE)
        self.play(Write(srect))
        self.wait()
        self.play(Write(rects))
        self.wait()
        for newrects in rect_list_cpy[1:]:
            self.play(Transform(rects, newrects))
        self.wait()
        self.play(FadeOut(rects))
        self.wait()
        graph_grp = VGroup()
        tgraph_grp = VGroup(
            self.get_graph(lambda t: t, x_min = 0, x_max = 2),
            self.get_graph(lambda t: t ** (1 / 2), x_min = 0, x_max = 2),
            self.get_graph(lambda t: np.log(1 + t), x_min = 0, x_max = 2),
            self.get_graph(lambda t: np.exp(t / 2) - 1, x_min = 0, x_max = 2),
            self.get_graph(lambda t: (9 - (t  - 2) ** 2) ** 0.5 - 5 ** 0.5, x_min = 0, x_max = 2),
        )
        ttgraph_grp = VGroup(
            self.get_graph(lambda t: 1 + t, x_min = 0, x_max = 2),
            self.get_graph(lambda t: 1 + t ** (1 / 2), x_min = 0, x_max = 2),
            self.get_graph(lambda t: 1 + np.log(1 + t), x_min = 0, x_max = 2),
            self.get_graph(lambda t: 1 + np.exp(t / 2) - 1, x_min = 0, x_max = 2),
            self.get_graph(lambda t: 1 + (9 - (t  - 2) ** 2) ** 0.5 - 5 ** 0.5, x_min = 0, x_max = 2),
        )
        cava_grp = VGroup()
        b_line = Line(15 * LEFT, 15 * RIGHT).set_color(GREY)
        self.play(Write(b_line))
        for n in range(1, 6):
            t = tgraph_grp[n - 1].copy()
            tt = ttgraph_grp[n - 1].copy()
            t = xyreflect(t)
            tt = xyreflect(tt)
            graph_grp.add(t.shift(5 * RIGHT).set_color(YELLOW))
            temp = cava_rectanges(t, tt, x_max = 2, stroke_width = 1 / 2 ** n, dx = 1 / 2 ** n, ashift = 0)
            temp.shift(4 * RIGHT)
            cava_grp.add(temp)
        rect_list = cava_list(base_gr, hei_gr, ashift = 1, x_max = 2, n_iterations = 5)
        rect_list_cpy = rect_list.copy()
        rects = rect_list_cpy[0]
        self.play(Write(rects))
        self.play(Write(graph_grp[0]), ReplacementTransform(rects.copy(), cava_grp[0]))
        self.play(FadeOut(graph_grp[0]), FadeOut(cava_grp[0]))
        for k, g, r in zip(graph_grp[1:], cava_grp[1:], rect_list_cpy[1:]):
            self.play(Transform(rects, r))
            self.wait()
            self.play(Write(k), ReplacementTransform(rects.copy(), g))
            self.wait()
            self.play(FadeOut(k), FadeOut(g))
            self.wait()
        arshape = VGroup(
            Line(self.input_to_graph_point(0, tgraph_grp[-1]), self.input_to_graph_point(0, ttgraph_grp[-1])),
            tgraph_grp[-1],
            Line(self.input_to_graph_point(2, tgraph_grp[-1]), self.input_to_graph_point(2, ttgraph_grp[-1])),
            ttgraph_grp[-1]).set_color(BLUE)
        arshape = xyreflect(arshape).shift(4 * RIGHT)
        lline = Line(self.input_to_graph_point(h.get_value(), tgraph_grp[-1]), self.input_to_graph_point(h.get_value(), ttgraph_grp[-1]))
        lline = xyreflect(lline).shift(4 * RIGHT)
        oline = Line(h.get_value() * UP, RIGHT + h.get_value() * UP)
        line_gp = VGroup(lline, oline).set_color(RED)
        self.play(Write(line_gp))
        self.wait()
        def lupd(obj):
            new_obj = Line(self.input_to_graph_point(h.get_value(), tgraph_grp[-1]), self.input_to_graph_point(h.get_value(), ttgraph_grp[-1]))
            new_obj = xyreflect(new_obj).shift(4 * RIGHT)
            n_obj = Line(h.get_value() * UP, RIGHT + h.get_value() * UP)
            n_gp = VGroup(new_obj, n_obj).set_color(RED)
            obj.become(n_gp)
        line_gp.add_updater(lupd)
        self.add(line_gp)
        self.play(Write(arshape))
        self.play(h.increment_value, 2, run_time = 4)
        self.play(h.increment_value, -2, run_time = 4)
        line_gp.clear_updaters()
        self.play(*[FadeOut(i) for i in [line_gp, rects, arshape, srect]], self.camera_frame.shift, RIGHT, titl.shift, RIGHT)
        gr1 = self.get_graph(lambda t: t ** 2, x_min = 0, x_max = 2)
        gr2 = self.get_graph(lambda t: t ** (1 / 2), x_min = 0, x_max = 2)
        graph_grp = VGroup(gr1, gr2.shift(UP))
        graph_grp = xyreflect(graph_grp)
        self.play(Write(graph_grp))
        rect_list = cava_list(gr1, gr2, ashift = 1, n_iterations = 5, x_max = 2)
        rects = rect_list[0]
        self.play(FadeIn(rects, run_time = 2))
        for new_rects in rect_list[1:]:
            self.transform_between_riemann_rects(rects, new_rects)
        gr3 = self.get_graph(lambda t: 0, x_min = 0, x_max = 2)
        gr4 = self.get_graph(lambda t: 1 + t ** 0.5 - t ** 2, x_min = 0, x_max = 2)
        ngr_grp = xyreflect(VGroup(gr3, gr4))
        ngr_grp.shift(6 * RIGHT)
        self.play(Write(ngr_grp))
        nrect_list = cava_list(gr3, gr4, ashift = 0, n_iterations = 5, x_max = 2)
        for n in nrect_list:
            n.shift(6 * RIGHT)
        nrects = nrect_list[0]
        self.play(FadeIn(nrects, run_time = 2))
        for nnew_rects in nrect_list[1:]:
            self.transform_between_riemann_rects(nrects, nnew_rects)
        lline = Line(self.input_to_graph_point(h.get_value(), gr1), self.input_to_graph_point(h.get_value(), gr2) + UP)
        lline = xyreflect(lline)
        oline = Line(self.input_to_graph_point(h.get_value(), gr3), self.input_to_graph_point(h.get_value(), gr4))
        oline = xyreflect(oline).shift(6 * RIGHT)
        line_gp = VGroup(lline, oline).set_color(RED)
        def llupd(obj):
            new_obj = Line(self.input_to_graph_point(h.get_value(), gr1), self.input_to_graph_point(h.get_value(), gr2) + UP)
            new_obj = xyreflect(new_obj)
            n_obj = Line(self.input_to_graph_point(h.get_value(), gr3), self.input_to_graph_point(h.get_value(), gr4))
            n_obj = xyreflect(n_obj).shift(6 * RIGHT)
            n_gp = VGroup(new_obj, n_obj).set_color(RED)
            obj.become(n_gp)
        self.play(Write(line_gp))
        line_gp.add_updater(llupd)
        self.add(line_gp)
        self.play(h.increment_value, 2, run_time = 4)
        self.play(h.increment_value, -2, run_time = 4)
        line_gp.clear_updaters()
        self.play(*[FadeOut(i) for i in [line_gp, rects, nrects, graph_grp, ngr_grp, titl, b_line]])
        self.wait(2)

class AreaByCavalieri(GraphScene, MovingCameraScene):
    CONFIG = {
        "x_min": -10,
        "x_max": 10,
        "x_axis_width": 20,
        "x_tick_frequency": 1,
        "x_leftmost_tick": None,  # Change if different from x_min
        "x_labeled_nums": None,
        "x_axis_label": "",
        "y_min": -5,
        "y_max": 5,
        "y_axis_height": 10,
        "y_tick_frequency": 1,
        "y_bottom_tick": None,  # Change if different from y_min
        "y_labeled_nums": None,
        "y_axis_label": "",
        "axes_color": GREY,
        "graph_origin": ORIGIN,
        "exclude_zero_label": True,
        "default_graph_colors": [BLUE, GREEN, YELLOW],
        "default_derivative_color": GREEN,
        "default_input_color": YELLOW,
        "default_riemann_start_color": BLUE,
        "default_riemann_end_color": GREEN,
        "area_opacity": 0.8,
        "num_rects": 50,
        "dot_kwargs": {
            "radius": 0.05,
        },
        "line_kwargs": {
            "stroke_width": 2,
        },
        "fill_triangle_kwargs": {
            # "fill_color": BLUE,
            "fill_opacity": .5,
            "stroke_width": 0,
        },
    }

    def setup(self):
        GraphScene.setup(self)
        MovingCameraScene.setup(self)
    
    def construct(self):
        self.setup_axes()
        self.play(FadeOut(self.axes))
        def cycloidpts(t, radius = 1, frac = 1, inverted = False):
            return (radius * t - (-1 if inverted else 1) * frac * radius * np.sin(t)) * RIGHT + (radius - radius * frac * np.cos(t)) * UP
        f = ValueTracker(1)
        theta = ValueTracker(0)
        b_line = Line(10 * LEFT, 10 * RIGHT).set_color(GREY)
        cyc = Cycloid(frac = f.get_value(), point_a = 4 * PI * LEFT, end_theta = 8 * PI)
        gen_grp = VGroup(Circle(), Circle(radius = f.get_value()).fade(0.5), Line(ORIGIN, f.get_value() * DOWN), Dot().shift(f.get_value() * DOWN)).set_color(RED)
        gen_grp.rotate(-theta.get_value()).move_to(theta.get_value() * RIGHT + UP)
        self.play(Write(b_line), Write(cyc), Write(gen_grp))
        archarea =  VMobject(fill_color = PURPLE, fill_opacity = 0.5, stroke_width = 0)
        archarea_pts = [cycloidpts(i, frac = f.get_value()) for i in np.arange(0, 2 * PI + 0.1, 0.1)]
        archarea.set_points_as_corners(archarea_pts)
        self.play(cyc.fade, 0.875, Write(archarea))
        self.wait()
        self.play(self.camera_frame.scale, 0.75, self.camera_frame.shift, PI * RIGHT + UP)
        titl = TexMobject("\\text{\\underline{Area under the Cycloidal Arch}}").move_to(PI * RIGHT + 3.5 * UP)
        cyc_grp = VGroup(cyc, archarea)
        def area_upd(obj):
            temp_c = Cycloid(frac = f.get_value(), point_a = 4 * PI * LEFT, end_theta = 8 * PI).fade(0.875)
            t_archarea =  VMobject(fill_color = PURPLE, fill_opacity = 0.5, stroke_width = 0)
            t_archarea_pts = [cycloidpts(i, frac = f.get_value()) for i in np.arange(0, 2 * PI + 0.1, 0.1)]
            t_archarea_pts = [ORIGIN] + t_archarea_pts + [2 * PI * RIGHT]
            t_archarea.set_points_as_corners(t_archarea_pts)
            t_grp = VGroup(temp_c, t_archarea)
            obj.become(t_grp)
        def gen_upd(obj):
            fr = f.get_value()
            t_grp = VGroup(Circle(), Circle(radius = f.get_value()).fade(0.5), Line(ORIGIN, fr * DOWN), Dot().shift(fr * DOWN)).set_color(RED)
            t_grp.rotate(-theta.get_value()).move_to(theta.get_value() * RIGHT + UP)
            obj.become(t_grp)
        cyc_grp.add_updater(area_upd)
        gen_grp.add_updater(gen_upd)
        self.add(cyc_grp)
        self.add(gen_grp)
        self.play(Write(titl))
        self.play(f.increment_value, -0.4, run_time = 3)
        one_arch_cyc = Cycloid(frac = f.get_value())
        self.play(theta.increment_value, 2 * PI, ShowCreation(one_arch_cyc), run_time = 5, rate_func = smooth)
        self.play(theta.increment_value, -2 * PI, run_time = 5, rate_func = smooth)
        cyc_grp.clear_updaters()
        harcharea = VMobject(fill_color = PURPLE, fill_opacity = 0.5, stroke_width = 0)
        harchareapts = [cycloidpts(i, frac = f.get_value()) for i in np.arange(0, PI, 0.1)]
        harchareapts += [cycloidpts(PI, frac = f.get_value())]
        harchareapts = [ORIGIN] + harchareapts + [PI * RIGHT]
        harcharea.set_points_as_corners(harchareapts)
        bbox = VMobject().set_color(YELLOW)
        bbox.set_points_as_corners([ORIGIN, PI * RIGHT, PI * RIGHT + 2 * UP, 2 * UP, ORIGIN])
        self.play(self.camera_frame.shift, PI / 2 * LEFT, titl.shift, PI / 2 * LEFT)
        self.play(Write(bbox))
        self.play(Transform(archarea, harcharea))
        self.wait()
        self.play(FadeOut(archarea))
        alp = ValueTracker(0)
        agen_grp = VGroup(Circle(), Circle(radius = f.get_value()).fade(1), Line(ORIGIN, f.get_value() * DOWN), Dot().shift(f.get_value() * DOWN)).set_color(RED)
        agen_grp.rotate(alp.get_value()).move_to(alp.get_value() * RIGHT + UP)
        t_line = Line(10 * LEFT + 2 * UP, 10 * RIGHT + 2 * UP).set_color(GRAY)
        self.play(Write(t_line))
        self.add(bbox)
        def agen_upd(obj):
            fr = f.get_value()
            t_grp = VGroup(Circle(), Circle(radius = f.get_value()).fade(1), Line(ORIGIN, fr * DOWN), Dot().shift(fr * DOWN)).set_color(RED)
            t_grp.rotate(alp.get_value()).move_to(alp.get_value() * RIGHT + UP)
            obj.become(t_grp)
        agen_grp.add_updater(agen_upd)
        self.add(agen_grp)
        acyc = Cycloid(frac = f.get_value(), inverted = True)
        self.play(alp.increment_value, 2 * PI, ShowCreation(acyc), run_time = 5)
        self.play(alp.increment_value, -2 * PI, run_time = 5)
        half_arch_area = VMobject(fill_color = BLUE, fill_opacity = 0.5, stroke_width = 0)
        lhalf_arch_area = VMobject(fill_color = BLUE, fill_opacity = 0.5, stroke_width = 0)
        half_arch_area_pts = [cycloidpts(i, frac = f.get_value()) for i in np.arange(0, 2 * PI + 0.1, 0.1)] + list(reversed([cycloidpts(i, frac = f.get_value(), inverted = True) for i in np.arange(0, 2 * PI + 0.1, 0.1)]))
        lhalf_arch_area_pts = [cycloidpts(i, frac = f.get_value()) for i in np.arange(0, PI + 0.1, 0.1)] + list(reversed([cycloidpts(i, frac = f.get_value(), inverted = True) for i in np.arange(0, PI + 0.1, 0.1)]))
        half_arch_area.set_points_as_corners(half_arch_area_pts)
        lhalf_arch_area.set_points_as_corners(lhalf_arch_area_pts)
        self.play(Write(half_arch_area))
        self.wait()
        self.play(FadeOut(half_arch_area), FadeIn(lhalf_arch_area))
        self.wait()
        small_circle = Circle(radius = f.get_value()).set_color(RED).shift(UP)
        cline = Line(gen_grp[-1].get_center(), agen_grp[-1].get_center()).set_color(GREEN)
        lcline = cline.copy().shift(PI * LEFT)
        line_grp = VGroup(cline, lcline)
        self.play(Write(line_grp))
        def lupd(obj):
            new_cl = Line(gen_grp[-1].get_center(), agen_grp[-1].get_center()).set_color(GREEN)
            new_lcl = new_cl.copy().shift(PI * LEFT + theta.get_value() * LEFT)
            new_obj = VGroup(new_cl, new_lcl)
            obj.become(new_obj)
        line_grp.add_updater(lupd)
        self.add(line_grp)
        self.play(
            ReplacementTransform(small_circle.copy(), small_circle.shift(PI * LEFT)),
            self.camera_frame.shift, PI / 2 * LEFT,
            titl.shift, PI / 2 * LEFT)
        self.play(b_line.fade, 0.75, t_line.fade, 0.75)
        for _ in range(4):
            self.play(theta.increment_value, PI / 4, alp.increment_value, PI / 4, run_time = 3)
        self.play(theta.increment_value, -PI, alp.increment_value, -PI, run_time = 6)
        self.wait()
        gen_grp.clear_updaters()
        agen_grp.clear_updaters()
        cline.clear_updaters()
        a_brace = Brace(Line(UP, 2 * UP), direction = LEFT, buff = 0.05, width_multiplier = 2.5)
        a_label = TexMobject("a").scale(0.5).next_to(a_brace, direction = LEFT, buff = 0.05)
        b_brace = Brace(Line(UP, UP + f.get_value() * DOWN), direction = LEFT, buff = 0.05, width_multiplier = 2.5)
        b_label = TexMobject("b").scale(0.5).next_to(b_brace, direction = LEFT, buff = 0.05)
        len_brace = Brace(Line(2 * UP, 2 * UP + PI * RIGHT), direction = UP, buff = 0.05, width_multiplier = 2.5)
        len_label = TexMobject("\\pi a").scale(0.5).next_to(len_brace, direction = UP, buff = 0.05)
        bre_brace = Brace(Line(2 * UP + PI * RIGHT, PI * RIGHT), direction = RIGHT, buff = 0.05, width_multiplier = 2.5)
        bre_label = TexMobject("2a").scale(0.5).next_to(bre_brace, direction = RIGHT, buff = 0.05)
        self.play(*[Write(i) for i in [a_brace, a_label, b_brace, b_label]], self.camera_frame.shift, DOWN)
        self.wait()
        self.play(*[Write(i) for i in [len_brace, len_label, bre_brace, bre_label]])
        sbox = VMobject().set_color(ORANGE)
        sbox.set_points_as_corners([
            UP + f.get_value() * DOWN,
            UP + f.get_value() * DOWN + PI * RIGHT,
            UP + f.get_value() * UP + PI * RIGHT,
            UP + f.get_value() * UP,
            UP + f.get_value() * DOWN,])
        sc_haarea = lhalf_arch_area.copy().scale(0.5).move_to(DOWN + PI * LEFT + LEFT)
        self.play(Write(sbox))
        sbox_brace = Brace(Line(PI * RIGHT + UP + f.get_value() * DOWN, PI * RIGHT + UP + f.get_value() * UP), direction = RIGHT, buff = 0.05)
        sbox_label = TexMobject("2b").scale(0.5).next_to(sbox_brace, direction = RIGHT, buff = 0.05)
        self.play(bre_brace.shift, 0.75 * RIGHT, bre_label.shift, 0.75 * RIGHT, Write(sbox_brace), Write(sbox_label))
        self.wait()
        lower_arch_area = VMobject(fill_color = ORANGE, fill_opacity = 0.5, stroke_width = 0)
        lower_arch_area_pts = [cycloidpts(i, frac = f.get_value(), inverted = True) for i in np.arange(0, PI + 0.1, 0.1)] + [PI * RIGHT + UP + f.get_value() * DOWN]
        lower_arch_area_pts += [lower_arch_area_pts[0]]
        lower_arch_area.set_points_as_corners(lower_arch_area_pts)
        upper_arch_area = lower_arch_area.copy().rotate(PI, about_point = PI / 2 * RIGHT + UP)
        self.play(Write(lower_arch_area), Write(upper_arch_area))
        self.wait()
        self.play(ReplacementTransform(lhalf_arch_area.copy(), sc_haarea))
        pib2 = TexMobject("=\\pi b^2").scale(0.5).next_to(sc_haarea, buff = 0)
        self.play(Write(pib2))
        self.wait()
        sc_haarea_cpy = sc_haarea.copy().shift(PI * RIGHT)
        twotimes = TexMobject("+\\text{  }2\\text{ }\\cross").scale(0.5).next_to(sc_haarea_cpy, buff = 0)
        lwrcopy = lower_arch_area.copy().scale_in_place(0.5).next_to(twotimes, buff = 0)
        twopiab = TexMobject("=2\\pi a b").scale(0.5).next_to(lwrcopy, buff = 0.1)
        self.play(Write(VGroup(sc_haarea_cpy, twotimes, lwrcopy, twopiab)))
        self.wait()
        self.play(*[FadeOut(i) for i in [sc_haarea_cpy, twopiab, twotimes]])
        self.play(lwrcopy.shift, 3 * LEFT)
        lwrcopy_area = TexMobject("=b(2a-b)\\pi / 2").scale(0.5).next_to(lwrcopy, buff = 0.1)
        self.play(Write(lwrcopy_area))
        brec = VMobject(fill_color = YELLOW, fill_opacity = 0.5, stroke_width = 0)
        brec.set_points_as_corners([
            ORIGIN,
            PI * RIGHT,
            UP + f.get_value() * DOWN + PI * RIGHT,
            UP + f.get_value() * DOWN,
            ORIGIN,])
        u_brec = brec.copy().shift(UP + f.get_value() * UP)
        self.play(Write(brec), Write(u_brec))
        brec_cpy = brec.copy().scale(0.5).move_to(PI * RIGHT + LEFT + DOWN)
        brec_area = TexMobject("=\\pi a(a - b)").scale(0.5).next_to(brec_cpy, buff = 0.1)
        self.play(ReplacementTransform(brec.copy(), brec_cpy), Write(brec_area))
        l_cpy = harcharea.copy().scale(0.5).move_to(2.25 * DOWN + PI * LEFT)
        area_gp = VGroup(
            TexMobject("=").scale(0.5),
            sc_haarea.copy(),
            TexMobject("+").scale(0.5),
            lwrcopy.copy(),
            TexMobject("+").scale(0.5),
            brec_cpy.copy()
        ).arrange(direction = RIGHT, buff = 0.25)
        area_gp.next_to(l_cpy, buff = 0.25)
        self.play(Write(l_cpy), Write(area_gp))
        self.wait()
        self.play(
            l_cpy.shift, PI / 2 * RIGHT,
            Transform(area_gp, TexMobject("=\\frac{\\pi}{2}(2a^2+b^2)").scale(0.5).move_to(0.125 * RIGHT + 2.25 * DOWN))
        )
        self.wait()
        self.play(
            *[FadeOut(i) for i in [l_cpy, lwrcopy, area_gp, pib2, sc_haarea, brec_cpy, brec_area, lwrcopy_area, a_brace, a_label, b_brace, b_label, len_brace, len_label, bre_brace, bre_label, sbox_brace, sbox_label]],
            self.camera_frame.move_to, PI * RIGHT + UP,
            titl.shift, PI * RIGHT)
        archarea =  VMobject(fill_color = PURPLE, fill_opacity = 0.5, stroke_width = 0)
        archarea_pts = [cycloidpts(i, frac = f.get_value()) for i in np.arange(0, 2 * PI + 0.1, 0.1)]
        archarea_pts = [ORIGIN] + archarea_pts + [2 * PI * RIGHT]
        archarea.set_points_as_corners(archarea_pts)
        self.play(*[FadeOut(i) for i in [brec, u_brec, bbox, lower_arch_area, upper_arch_area, lhalf_arch_area, sbox, acyc]], FadeIn(archarea))
        self.wait()
        archarea_label = TexMobject("\\pi (2a^2+b^2)").shift(PI * RIGHT + DOWN)
        self.play(Write(archarea_label))
        self.wait(2)

if __name__ == "__main__":
    module_name = os.path.basename(__file__)
    #command_A = "manim -p  --video_dir ~/Downloads/  "
    clear_cmd = "cls"
    command_A = "python -m manim "
    command_B = module_name + " " + "LearningCava" + " -p -n 61"
    #command_B = module_name + " -p"
    os.system(clear_cmd)
    os.system(command_A + command_B)