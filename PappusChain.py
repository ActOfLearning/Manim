from manimlib.imports import *
import numpy as np

def Circle_Inversion(obj, about_point = ORIGIN, inversion_radius = 1.0):
    def inverse_of_a_point(pt0):
        norm_of_pt = np.dot(pt0, pt0)
        if norm_of_pt == 0:
            norm_of_pt += 0.0001
        return pt0 / norm_of_pt
    def inverse_of_a_circle(obj):
        center_of_the_circle = obj.get_center()
        radius_of_the_circle = np.linalg.norm(center_of_the_circle - obj.get_start())
        inversion_factor = np.dot(center_of_the_circle, center_of_the_circle) - radius_of_the_circle ** 2
        if inversion_factor == 0:
            res_object = Line(10 * UP, 10 * DOWN)
            res_object.move_to(RIGHT * (1 / 2 / radius_of_the_circle))
            rot_angle = np.angle(center_of_the_circle[0] + center_of_the_circle[1] * 1j)
            res_object.rotate_about_origin(rot_angle)
            return res_object
        inversion_factor = 1.0 / inversion_factor
        return Circle(radius = abs(inversion_factor) * radius_of_the_circle).move_to(inversion_factor * center_of_the_circle).set_color(TEAL)
    def circum_center(pt0):
        pt = pt0[:]
        side_a, side_b, side_c = pt[2] - pt[1], pt[0] - pt[2], pt[1] - pt[0]
        side_a, side_b, side_c = np.dot(side_a, side_a), np.dot(side_b, side_b), np.dot(side_c, side_c)
        circum_x, circum_y, circum_z = side_a * (side_b + side_c - side_a), side_b * (side_c + side_a - side_b), side_c * (side_a + side_b - side_c)
        return (circum_x * pt[0] + circum_y * pt[1] + circum_z * pt[2]) / (circum_x + circum_y + circum_z)
    def inverted_arc(c0, r0, ang0):
        return Arc(ang0[0], ang0[1] - ang0[0] + (2 * PI if ang0[1] < ang0[0] else 0), radius = r0, arc_center = c0)
    def arg_of_a_point(pt0):
        return np.angle(pt0[0] + 1j * pt0[1])
    def inverse_of_a_line(obj):
        line_st, line_en = obj.get_start(), obj.get_end()
        inv_line_start = inverse_of_a_point(line_st)
        inv_line_end = inverse_of_a_point(line_en)
        if abs(np.linalg.det([line_st, line_en, OUT])) >= 0.001:
            circum = circum_center([ORIGIN, inv_line_start, inv_line_end])
            circum_norm = np.linalg.norm(circum)
            temp_mat = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]])
            temp_mat = np.dot(np.dot(temp_mat, np.transpose(line_st)), line_en)
            return inverted_arc(circum, circum_norm, [arg_of_a_point(inv_line_start - circum), arg_of_a_point(inv_line_end - circum)] if temp_mat > 0 else [arg_of_a_point(inv_line_end - circum), arg_of_a_point(inv_line_start - circum)])
        else:
            if np.dot(line_st, line_en) > 0:
                return Line(inv_line_start, inv_line_end)
            else:
                return VGroup(Line(inv_line_start, 6 * inv_line_start / (0.1 + np.linalg.norm(inv_line_start))), Line(inv_line_end, 6 * inv_line_end / (0.1 + np.linalg.norm(inv_line_end))))
    def inverse_of_a_polygon(obj):
        polygon_vertices = obj.get_points()
        inverted_points = VGroup()
        for i in range(1, len(polygon_vertices)):
            inverted_points.add(inverse_of_a_line(Line(polygon_vertices[i], polygon_vertices[i - 1])))
        return inverted_points.set_color(TEAL)
    shifted_obj = obj.copy()
    shifted_obj.shift(-about_point)
    shifted_obj.scale(1 / inversion_radius, about_point = ORIGIN)
    if type(obj).__name__ == "Circle":
        res_obj = inverse_of_a_circle(shifted_obj)
    elif type(obj).__name__ == "Dot":
        res_obj = Dot(radius = DEFAULT_DOT_RADIUS / inversion_radius).move_to(inverse_of_a_point(shifted_obj.get_center()))
    elif type(obj).__name__ == "VGroup":
        res_obj = VGroup()
        for objl in shifted_obj:
            res_obj.add(Circle_Inversion(objl))
    else:
        res_obj = inverse_of_a_polygon(shifted_obj)
    res_obj.scale(inversion_radius, about_point = ORIGIN)
    res_obj.shift(about_point)
    return res_obj.set_color(TEAL)

def ParallalLine(pt0, line0):
    res = VGroup()
    if type(pt0).__name__ == "Dot" and type(line0).__name__ == "Line":
        line_start, line_end = line0.get_start_and_end()
        p_cent = pt0.get_center()
        if np.linalg.norm(line_end - line_start) > 0:
            normalized_line = (line_end - line_start) / np.linalg.norm(line_end - line_start)
            return Line(p_cent, p_cent + normalized_line)
    return res

def PerpendicularFoot(pt0, line0):
    res = VGroup()
    if type(pt0).__name__ == "Dot" and type(line0).__name__ == "Line":
        line_start, line_end = line0.get_start_and_end()
        p_cent = pt0.get_center()
        if np.linalg.norm(line_end - line_start) > 0:
            normalized_line = (line_end - line_start) / np.linalg.norm(line_end - line_start)
            line_start += np.dot(p_cent - line_start, normalized_line) * normalized_line
            return Dot().move_to(line_start)
    return res

def PerpendicularLine(pt0, line0):
    line_start, line_end = line0.get_start_and_end()
    line_unit = (line_end - line_start) / np.linalg.norm(line_end - line_start)
    perp_unit = np.cross(line_unit, OUT)
    perp_foot = PerpendicularFoot(pt0, line0).get_center()
    return Line(perp_foot, perp_unit + perp_foot)

def DottoCircle(obj1):
    cen1, rad1 = ORIGIN, 0.0
    if type(obj1).__name__ == "Dot" or type(obj1).__name__ == "Circle":
        cen1 = obj1.get_center()
        if type(obj1).__name__ == "Dot":
            rad1 = 0.0
        else:
            rad1 = np.linalg.norm(cen1 - obj1.get_start())
    return cen1, rad1

def RadicalAxis(obj1, obj2):
    if type(obj1).__name__ == "Dot" or type(obj1).__name__ == "Circle":
        cen1, rad1 = DottoCircle(obj1)
    if type(obj2).__name__ == "Dot" or type(obj2).__name__ == "Circle":
        cen2, rad2 = DottoCircle(obj2)
    dist_cent = np.linalg.norm(cen1 - cen2)
    if dist_cent > 0:
        res = rad1 ** 2 - rad2 ** 2
        res = dist_cent + res / dist_cent
        res /= 2
        res /= dist_cent
        res = cen1 + res * (cen2 - cen1)
        return PerpendicularLine(Dot().move_to(res), Line(cen1, cen2))

def Intersection(*args):
    res_obj = VGroup()
    obj1, obj2 = args[0], args[1]
    if type(obj1).__name__ == "Line":
        obj1_start, obj1_end = obj1.get_start_and_end()
        ab = obj1_end - obj1_start
        ab_unit = ab / np.linalg.norm(ab)
        if type(obj2).__name__ == "Line":
            obj2_start, obj2_end = obj2.get_start_and_end()
            cd = obj2_end - obj2_start
            ac = obj2_start - obj1_start
            abcd = np.cross(ab, cd)
            abcd_n = np.dot(abcd, abcd)
            if abcd_n != 0:
                abcd = abcd / abcd_n
                abcd = np.dot(np.cross(ac, cd), abcd)
                abcd = obj1_start + ab * abcd
                res_obj.add(Dot(abcd))
        elif type(obj2).__name__ == "Circle":
            cen = Dot().move_to(obj2.get_center())
            p_foot = PerpendicularFoot(cen, obj1)
            p_dist = np.linalg.norm(cen.get_center() - p_foot.get_center())
            rad = np.linalg.norm(obj2.get_center() - obj2.get_start())
            if p_dist == rad:
                res_obj.add(p_foot)
            elif p_dist < rad:
                chord_dist = rad * rad - p_dist * p_dist
                chord_dist = np.sqrt(chord_dist)
                res_obj.add(Dot().move_to(p_foot.get_center() + chord_dist * ab_unit), Dot().move_to(p_foot.get_center() - chord_dist * ab_unit))
    elif type(obj1).__name__ == "Circle":
        cen1, cen2 = obj1.get_center(), obj2.get_center()
        dist_cen = np.linalg.norm(cen1 - cen2)
        rad1, rad2 = np.linalg.norm(cen1 - obj1.get_start()), np.linalg.norm(cen2 - obj2.get_start())
        if abs(dist_cen - rad1 - rad2) < 1e-3:
            res_obj.add(Dot().move_to((rad2 * cen1 + rad1 * cen2) / (rad1 + rad2)))
        elif dist_cen < rad1 + rad2:
            p_foot = RadicalAxis(obj1, obj2)
            #p_foot = PerpendicularLine(p_foot, Line(cen1, cen2))
            return Intersection(p_foot, obj1)
    return res_obj
            
def Midpoint(*args):
    if type(args[0]).__name__ == "Line":
        line_start, line_end = args[0].get_start_and_end()
        mid_pt = (line_start + line_end) / 2
        mid_pt = Dot().move_to(mid_pt)
        return mid_pt
    elif type(args[0]).__name__ == "Dot" and type(args[1]).__name__ == "Dot":
        mid_pt = (args[0].get_center() + args[1].get_center()) / 2
        mid_pt = Dot().move_to(mid_pt)
        return mid_pt

def PerpendicularBisector(*args):
    mid_pt = Midpoint(*args)
    if type(args[0]).__name__ == "Line":
        print("It's a line..")
        return PerpendicularLine(mid_pt, args[0])
    elif type(args[0]).__name__ == "Dot" and type(args[1]).__name__ == "Dot":
        connec_line = Line(args[0].get_center(), args[1].get_center())
        return PerpendicularLine(mid_pt, connec_line)

def AngularBisector(*args):
    if type(args[0]).__name__ == "Line" and type(args[1]).__name__ == "Line":
        intersec = Intersection(*args)
        if len(list(intersec)) > 0:
            return AngularBisector(Dot().move_to(args[0].get_end()), intersec, Dot().move_to(args[1].get_end()))
        else:
            return VGroup()
    elif type(args[0]).__name__ == "Dot" and type(args[1]).__name__ == "Dot" and type(args[2]).__name__ == "Dot":
        pt_a = args[0].get_center()
        pt_b = args[1].get_center()
        pt_c = args[2].get_center()
        vec_ba = pt_a - pt_b
        vec_bc = pt_c - pt_b
        return Line(pt_b, vec_ba + vec_bc)

def RadicalCenter(*args):
    if type(args[0]).__name__ == "Circle" and type(args[1]).__name__ == "Circle" and type(args[2]).__name__ == "Circle":
        radic12 = RadicalAxis(args[0], args[1])
        radic23 = RadicalAxis(args[1], args[2])
        return Intersection(radic12, radic23)
    return Dot()

def HomotheticCenter(obj1, obj2):
    cen1, cen2 = obj1.get_center(), obj2.get_center()
    rad1, rad2 = np.linalg.norm(cen1 - obj1.get_start()), np.linalg.norm(cen2 - obj2.get_start())
    res_obj = VGroup()
    res_obj.add(Dot().move_to((rad1 * cen2 + rad2 * cen1) / (rad1 + rad2)))
    if rad1 != rad2:
        res_obj.add(Dot().move_to((rad1 * cen2 - rad2 * cen1) / (rad1 - rad2)))
    return res_obj

def CircleTangentLines(obj, pt):
    cen = obj.get_center()
    rad = np.linalg.norm(cen - obj.get_start())
    dist = np.linalg.norm(cen - pt.get_center())
    res_obj = VGroup()
    if dist < rad:
        return res_obj.add(pt)
    temp_c = Circle(arc_center = (cen + pt.get_center()) / 2, radius = dist / 2)
    temp_c = Intersection(temp_c, obj)
    for i in temp_c:
        res_obj.add(Line(pt.get_center(), i.get_center()))
    return res_obj

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
        self.setup_axes(animate = False)
        func = lambda x: x ** 3 + 1e-8
        parabola = self.get_graph(
            func,
            x_min = -50,
            x_max = 50,
            color = BLUE,
        )
        c = Circle()
        def Circle_Inverse_Func(funcc, x_min = -5, x_max = 5, step_size = 0.01, about_point = ORIGIN, inversion_radius = 1):
            res_vm = VMobject()
            res_list = [Circle_Inversion(Dot().move_to(np.array([x, funcc(x), 0])), about_point, inversion_radius).get_center() for x in np.arange(x_min, x_max, step_size)]
            res_vm.set_points_as_corners(res_list)
            return res_vm
        parabola_inv = Circle_Inverse_Func(func)
        self.play(ShowCreation(parabola), Write(c), ShowCreation(parabola_inv))
        self.wait(2)

class IntroductionScene(GraphScene, MovingCameraScene):
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
            "radius": 0.25,
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
        self.setup_axes(animate = False)
        tots = 2
        n = 2
        r = n * n / (2 * n * n - 1) if n == tots else (2 * n * n - 1 - (1 + 4 * n * n * (tots * tots - 1)) ** 0.5) / 2 / (n * n - tots * tots)
        rad1 = r * tots / 2
        rad2 = (1 - r) * tots / 2
        trad = rad1 + rad2
        c1, c2 = Circle(radius = rad1).move_to(rad1 * RIGHT).set_color(GREEN), Circle(radius = rad2).move_to(2 * rad1 * RIGHT + rad2 * RIGHT).set_color(LIGHT_PINK)
        c3 = Circle(radius = trad).move_to(trad * RIGHT).set_color(GREEN)
        self.camera_frame.save_state()
        self.play(self.camera_frame.scale, 0.375, self.camera_frame.move_to, RIGHT)
        r = rad1 / trad
        c_uplist, c_downlist = [], []
        center_list = []
        for m in range(- n * n - 2 * n + 2 - 5, 1 + n + 2 * n + 5):
            if m == 0:
                continue
            deno = 2 * (m * m * (1 - r) * (1 - r) + r)
            xm, ym = r * (1 + r) / deno, 2 * m * r * (1 - r) / deno
            rm = (1 - r) * r / deno
            xm, ym , rm = 2 * trad * xm, 2 * trad * ym, 2 * trad * rm
            if m < 0:
                c_downlist += [Circle(radius = rm).move_to(xm * RIGHT + ym * UP).set_color(LIGHT_PINK)]
            else:
                c_uplist += [Circle(radius = rm).move_to(xm * RIGHT + ym * UP).set_color(LIGHT_PINK)]
                center_list += [Dot(xm * RIGHT + ym * UP).scale(0.1).set_color(YELLOW)]
        def placelabel(label, obj, alpha = 0 * DEGREES, sc = 1):
            temp_dt = Dot(obj.get_start()).rotate(alpha, about_point = obj.get_center())
            label.scale(sc).next_to(temp_dt, direction = RIGHT, buff = 0)
        label_c1 = TextMobject("$C_A$")
        cen_c1 = Dot(c1.get_center()).scale(0.1)
        placelabel(label_c1, c1, alpha = -45 * DEGREES, sc = 0.2)
        label_c3 = TextMobject("$C_B$")
        cen_c3 = Dot(c3.get_center()).scale(0.1)
        placelabel(label_c3, c3, alpha = -45 * DEGREES, sc = 0.2)
        label_p0 = TextMobject("$P_0$").scale(0.2).next_to(c2.get_center(), direction = RIGHT, buff = 0)
        cen_p0 = Dot(c2.get_center()).scale(0.1)
        labellist = [TextMobject("$P_1$"), TextMobject("$P_2$"), TextMobject("$P_3$")]
        for c in range(3):
            labellist[c].scale(0.2).move_to(c_uplist[c].get_center())
        base_line = Line(c3.get_left() + LEFT / 2, c3.get_right() + RIGHT / 2).set_color(ORANGE)
        title_card = TexMobject("\\underline{\\text{Pappus Chain}}").scale(0.5).move_to(cen_c3.get_center() + 1.5 * UP)
        self.play(self.camera_frame.shift, 0.25 * UP, Write(title_card))
        self.wait(2)
        self.play(Write(c1), Write(c3))
        self.wait(2)
        self.play(Write(label_c1), Write(label_c3))
        self.wait(2)
        self.play(Write(c2))
        self.wait(2)
        self.play(Write(label_p0))
        self.wait(2)
        for c in range(3):
            self.play(Write(labellist[c]), Write(c_uplist[c]))
            self.wait(2)
        rem_vg = VGroup(*center_list[3:], *c_uplist[3:])
        self.play(Write(rem_vg))
        self.wait(2)
        self.play(Write(VGroup(*c_downlist)))
        self.wait(2)
        self.play(*[FadeOut(i) for i in labellist + [label_p0, label_c1, label_c3]])
        self.play(*[FadeIn(i) for i in center_list[:3] + [cen_c1, cen_p0]])
        self.wait(2)
        self.play(Write(base_line))
        self.wait(2)
        focus_circ = 0
        while focus_circ < 6:
            temp_c = c_uplist[focus_circ].copy()
            height_line = Line(temp_c.get_center(), PerpendicularFoot(Dot(temp_c.get_center()), base_line).get_center()).set_color(PURPLE)
            temp_vg = VGroup()
            temp_rad = np.linalg.norm(temp_c.get_center() - temp_c.get_start())
            temp_cnt, temp_dist = 0, temp_rad
            self.play(Write(height_line))
            while temp_cnt <= focus_circ:
                temp_vg.add(temp_c.copy().shift(temp_dist * DOWN))
                temp_dist += 2 * temp_rad
                temp_cnt += 1
            self.play(Transform(temp_c, temp_vg))
            self.wait(2)
            self.play(Transform(temp_c, c_uplist[focus_circ].copy()), FadeOut(height_line))
            self.remove(temp_c)
            focus_circ += 1
        self.wait(2)
        invcirt = Circle(radius = 1).move_to(4.5 * RIGHT).fade(1)
        cirt_grp = VGroup(c1, c2, c3, *c_uplist, *c_downlist)
        testc1_inv = Circle_Inversion(cirt_grp, invcirt.get_center()).fade(0.5)
        self.play(Write(testc1_inv))
        def testup(obj):
            temp_grp = VGroup(c1, c2, c3, *c_uplist, *c_downlist)
            obj.become(Circle_Inversion(temp_grp, invcirt.get_center()).fade(0.5))
        testc1_inv.add_updater(testup)
        self.add(testc1_inv)
        self.play(invcirt.move_to, 2.5 * LEFT, run_time = 10, rate_func = linear)
        testc1_inv.clear_updaters()
        self.play(FadeOut(testc1_inv), FadeOut(invcirt))
        self.wait(2)
        focus_circ = 2
        orig = Dot().scale(0.25)
        self.play(Write(orig))
        self.play(FadeOut(VGroup(*c_downlist)))
        t_lines = CircleTangentLines(c_uplist[focus_circ], orig)
        t_pt = t_lines[0].get_end()
        t_line = DashedLine(ORIGIN, 1.1 * t_pt).set_color(GRAY)
        self.play(Write(t_line))
        self.wait(2)
        ortho_arc = Arc(radius = t_lines[0].get_length(), start_angle = 80 * DEGREES, angle = -160 * DEGREES).set_color(GOLD).fade(0.5)
        self.play(Write(ortho_arc))
        self.wait(2)
        self.play(FadeOut(t_line))
        self.wait(2)
        c1_inv, c2_inv, c3_inv = Circle_Inversion(c1, inversion_radius = t_lines[0].get_length()), Circle_Inversion(c2, inversion_radius = t_lines[0].get_length()), Circle_Inversion(c3, inversion_radius = t_lines[0].get_length())
        c_uplist_inv = []
        c1_inter = Intersection(c1, Circle(radius = t_lines[0].get_length())).set_color(YELLOW)
        c3_inter = Intersection(c3, Circle(radius = t_lines[0].get_length())).set_color(YELLOW)
        for i in c1_inter:
            i.scale(0.1875)
        for i in c3_inter:
            i.scale(0.1875)
        for i in c_uplist:
            c_uplist_inv += [Circle_Inversion(i, inversion_radius = t_lines[0].get_length())]
        self.play(FadeInFromLarge(c1_inter, scale_factor = 3))
        self.wait(2)
        self.play(Write(c1_inv))
        self.wait(2)
        self.play(FadeOut(c1_inter))
        self.wait(2)
        self.play(FadeInFromLarge(c3_inter, scale_factor = 3))
        self.wait(2)
        self.play(Write(c3_inv))
        self.wait(2)
        self.play(FadeOut(c3_inter))
        self.wait(2)
        self.camera_frame.save_state()
        self.play(self.camera_frame.scale, 0.5, self.camera_frame.shift, 0.25 * UP)
        self.wait(2)
        ss1 = Dot().scale(0.1875).move_to(c1.get_start()).set_color(YELLOW)
        ss3 = Dot().scale(0.1875).move_to(c3.get_start()).set_color(YELLOW)
        s1 = Dot().scale(0.1875).move_to(Intersection(base_line, c1_inv).get_center()).set_color(YELLOW)
        s3 = Dot().scale(0.1875).move_to(Intersection(base_line, c3_inv).get_center()).set_color(YELLOW)
        self.play(FadeInFromLarge(ss1, scale_factor = 5))
        self.play(FadeInFromLarge(s1, scale_factor = 5))
        self.wait(2)
        self.play(FadeInFromLarge(ss3, scale_factor = 5))
        self.play(FadeInFromLarge(s3, scale_factor = 5))
        self.wait(2)
        inv_grp = VGroup(c2_inv)
        inv_label_ps = VGroup(
            TexMobject("P_0'").scale(0.15625).move_to(c2_inv.get_center()),
            TexMobject("P_1'").scale(0.15625).move_to(c_uplist_inv[0].get_center()),
            TexMobject("P_2'").scale(0.15625).move_to(c_uplist_inv[1].get_center()),
            TexMobject("P_3'").scale(0.15625).move_to(c_uplist_inv[2].get_center())
        )
        self.play(Write(c2_inv), Write(inv_label_ps[0]))
        self.wait(2)
        self.play(*[FadeOut(i) for i in [s1, s3, ss1, ss3]])
        self.wait(2)
        cnt = 1
        for i in c_uplist_inv[:focus_circ]:
            self.play(Write(i), Write(inv_label_ps[cnt]))
            self.wait(2)
            inv_grp.add(i)
            cnt += 1
        self.play(FadeOut(inv_label_ps))
        self.wait(2)
        height_line = Line(c_uplist[focus_circ].get_center(), PerpendicularFoot(Dot(c_uplist[focus_circ].get_center()), base_line).get_center()).set_color(PURPLE)
        foc_rad = np.linalg.norm(c_uplist_inv[0].get_center() - c_uplist_inv[0].get_start())
        self.play(inv_grp.shift, foc_rad * UP, FadeOut(ortho_arc))
        self.wait(2)
        self.play(Write(height_line))
        #inv_grp.shift(foc_rad * DOWN)
        self.wait(2)

class InversionIntroduction(GraphScene, MovingCameraScene):
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
        "default_graph_colors": [YELLOW, GREEN, YELLOW],
        "default_derivative_color": GREEN,
        "default_input_color": YELLOW,
        "default_riemann_start_color": YELLOW,
        "default_riemann_end_color": GREEN,
        "area_opacity": 0.8,
        "num_rects": 50,
        "dot_kwargs": {
            "radius": 0.25,
        },
        "line_kwargs": {
            "stroke_width": 2,
        },
        "fill_triangle_kwargs": {
            # "fill_color": YELLOW,
            "fill_opacity": .5,
            "stroke_width": 0,
        },
    }

    def setup(self):
        GraphScene.setup(self)
        MovingCameraScene.setup(self)
    
    def construct(self):
        self.setup_axes(animate = False)
        rad_of_inv = 2
        Inv_Circle = Circle(radius = rad_of_inv).set_color(PURPLE)
        Orig = Dot()
        Inv_center = TexMobject("O").scale(0.75)
        Inv_center.next_to(Orig, direction = DOWN + RIGHT, buff = 0.00)
        inv_label = VGroup(TexMobject("OA \\cdot OA'=r^2"), TexMobject("=OB \\cdot OB'")).scale(0.75)
        inv_label.arrange_submobjects(direction = RIGHT, buff = 0.1).shift(3.5 * UP)
        Apos = 0.9 * RIGHT
        A = TexMobject("A").scale(0.75).set_color(RED)
        Adot = Dot().move_to(Apos).set_color(RED)
        A.next_to(Adot, direction = DOWN + RIGHT, buff = 0.00)
        Adashdot = Circle_Inversion(Adot, inversion_radius = rad_of_inv)
        Adash = TexMobject("A'").scale(0.75).set_color(RED)
        Adash.next_to(Adashdot, direction = DOWN + RIGHT, buff = 0.00)
        rayOA = Line(ORIGIN, 10 * RIGHT).set_color(GREY)
        Bpos = 1 * UP + 1 * RIGHT
        B = TexMobject("B").scale(0.75).set_color(GREEN)
        Bdot = Dot().move_to(Bpos).set_color(GREEN)
        B.next_to(Bdot, direction = UL, buff = 0.00)
        Bdashdot = Circle_Inversion(Bdot, inversion_radius = rad_of_inv)
        Bdash = TexMobject("B'").scale(0.75).set_color(GREEN)
        Bdash.next_to(Bdashdot, direction = UL, buff = 0.00)
        rayOB = Line(ORIGIN, 10 * Bdot.get_center()).set_color(GREY)
        A_vg = VGroup(rayOA, A, Adot, Adash, Adashdot)
        def ptA_updater(obj):
            ar, a, ado, ad, addo = obj
            a.next_to(ado, direction = DOWN + RIGHT, buff = 0.00)
            new_addo = Circle_Inversion(ado, inversion_radius = rad_of_inv)
            addo.become(new_addo)
            ad.next_to(addo, direction = DOWN + RIGHT, buff = 0.00)
            new_ar = Line(ORIGIN, 10 * ado.get_center()).set_color(GREY)
            ar.become(new_ar)
        self.play(Write(Inv_Circle), Write(Inv_center), Write(Orig))
        self.wait(2)
        self.play(Write(Adot), Write(A))
        self.wait(2)
        self.play(Write(rayOA), Write(inv_label[0]))
        self.wait(2)
        self.play(Write(Adash), Write(Adashdot))
        self.wait(2)
        A_vg.add_updater(ptA_updater)
        self.add(A_vg)
        self.play(Adot.shift, 1.5 * RIGHT, rate_func = there_and_back, run_time = 3)
        self.wait(2)
        self.play(*[Write(i) for i in [B, Bdot, Bdash, Bdashdot, rayOB, inv_label[1]]])
        self.wait(2)
        A_vg.clear_updaters()
        AB = DashedLine(Adot.get_center(), Bdot.get_center())
        AdBd = DashedLine(Adashdot.get_center(), Bdashdot.get_center())
        self.play(Write(AB), Write(AdBd))
        self.wait(2)
        angleA = Sector(arc_center = Adot.get_center(), outer_radius = 0.375, start_angle = AB.get_angle(), angle = PI - AB.get_angle() + rayOA.get_angle(), fill_color = RED)
        angleB = Sector(arc_center = Bdot.get_center(), outer_radius = 0.375, start_angle = PI + rayOB.get_angle(), angle = AB.get_angle() - rayOB.get_angle(), fill_color = GREEN)
        angleAdash = Sector(arc_center = Adashdot.get_center(), outer_radius = 0.375, start_angle = AdBd.get_angle(), angle = PI - AdBd.get_angle() + rayOA.get_angle(), fill_color = GREEN)
        angleBdash = Sector(arc_center = Bdashdot.get_center(), outer_radius = 0.375, start_angle = PI + rayOB.get_angle(), angle = AdBd.get_angle() - rayOB.get_angle(), fill_color = RED)
        tri_vg = VGroup(rayOA, rayOB, AB, AdBd, A, Adot, Adash, Adashdot, B, Bdot, Bdash, Bdashdot, angleA, angleB, angleAdash, angleBdash)
        label_ratio_of_sides = TexMobject("\\frac{OA}{OB}=\\frac{OB'}{OA'}").scale(0.75)
        common_o = TexMobject("\\text{ and }\\angle O\\text{ is common}").scale(0.75)
        sim_triangles = TexMobject("\\implies \\triangle{OAB} \\equiv \\triangle{OB'A'}").scale(0.75).set_color(YELLOW)
        text_vg = VGroup(label_ratio_of_sides, common_o).arrange_submobjects(direction = RIGHT, buff = 0.2).shift(5 * LEFT + 2 * UP)
        self.play(ReplacementTransform(inv_label.copy(), text_vg[0]), self.camera_frame.shift, 1.5 * LEFT, inv_label.shift, 1.5 * LEFT)
        self.wait(2)
        self.play(Write(text_vg[1]))
        self.wait(2)
        sim_triangles.move_to(text_vg).shift(1.0 * DOWN)
        sim_triangles.align_to(text_vg[0], LEFT)
        self.play(Write(sim_triangles))
        self.wait(2)
        imp_similarity = VGroup(TexMobject("\\therefore \\angle A = \\angle B' \\text{ and } \\angle B = \\angle A',").scale(0.75),
        TexMobject("\\frac{B'A'}{AB}=\\frac{OB'}{OA}").scale(0.75))
        imp_similarity.arrange_submobjects(direction = DOWN, buff = 0.5)
        imp_similarity.align_to(text_vg[0], LEFT).shift(1.5 * DOWN)
        self.play(Write(imp_similarity[0]))
        self.wait(2)
        self.play(Write(angleA), Write(angleBdash))
        self.wait(2)
        self.play(Write(angleB), Write(angleAdash))
        self.wait(2)
        self.play(Write(imp_similarity[1]))
        self.wait(2)
        self.play(Transform(imp_similarity[1], TexMobject("A'B'=AB\\cdot\\frac{r^2}{OA \\cdot OB}").scale(0.75).move_to(imp_similarity[1].get_center())))
        self.wait(2)
        def group_updater(obj):
            ar, br, ab, adbd, a, ado, ad, addo, b, bdo, bd, bddo, anga, angb, angad, angbd = obj
            a.next_to(ado, direction = DOWN + RIGHT, buff = 0.00)
            new_addo = Circle_Inversion(ado, inversion_radius = rad_of_inv)
            addo.become(new_addo)
            ad.next_to(addo, direction = DOWN + RIGHT, buff = 0.00)
            new_ar = Line(ORIGIN, 10 * ado.get_center()).set_color(GREY)
            ar.become(new_ar)
            b.next_to(bdo, direction = UL, buff = 0.00)
            new_bddo = Circle_Inversion(bdo, inversion_radius = rad_of_inv)
            bddo.become(new_bddo)
            bd.next_to(bddo, direction = UL, buff = 0.00)
            new_br = Line(ORIGIN, 10 * bdo.get_center()).set_color(GREY)
            br.become(new_br)
            new_ab = DashedLine(ado.get_center(), bdo.get_center())
            new_adbd = DashedLine(addo.get_center(), bddo.get_center())
            ab.become(new_ab)
            adbd.become(new_adbd)
            new_anga = Sector(arc_center = ado.get_center(), outer_radius = 0.375, start_angle = ab.get_angle(), angle = PI - ab.get_angle() + ar.get_angle(), fill_color = RED)
            new_angb = Sector(arc_center = bdo.get_center(), outer_radius = 0.375, start_angle = PI + br.get_angle(), angle = ab.get_angle() - br.get_angle(), fill_color = GREEN)
            new_angad = Sector(arc_center = addo.get_center(), outer_radius = 0.375, start_angle = adbd.get_angle(), angle = PI - adbd.get_angle() + ar.get_angle(), fill_color = GREEN)
            new_angbd = Sector(arc_center = bddo.get_center(), outer_radius = 0.375, start_angle = PI + br.get_angle(), angle = adbd.get_angle() - br.get_angle(), fill_color = RED)
            anga.become(new_anga)
            angb.become(new_angb)
            angad.become(new_angad)
            angbd.become(new_angbd)
        tri_vg.add_updater(group_updater)
        self.add(tri_vg)
        self.play(Adot.shift, DOWN, Bdot.shift, 0.5 * RIGHT, run_time = 5, rate_func = there_and_back)
        self.wait(2)
        tri_vg.clear_updaters()
        self.wait(2)

class InversionCurves(GraphScene, MovingCameraScene):
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
        "default_graph_colors": [YELLOW, GREEN, YELLOW],
        "default_derivative_color": GREEN,
        "default_input_color": YELLOW,
        "default_riemann_start_color": YELLOW,
        "default_riemann_end_color": GREEN,
        "area_opacity": 0.8,
        "num_rects": 50,
        "dot_kwargs": {
            "radius": 0.25,
        },
        "line_kwargs": {
            "stroke_width": 2,
        },
        "fill_triangle_kwargs": {
            # "fill_color": YELLOW,
            "fill_opacity": .5,
            "stroke_width": 0,
        },
    }

    def setup(self):
        GraphScene.setup(self)
        MovingCameraScene.setup(self)
    
    def construct(self):
        self.setup_axes(animate = False)
        def Circle_Inverse_Func(grap, x_min = -10, x_max = 10, step_size = 0.01, about_point = ORIGIN, inversion_radius = 1):
            res_vm = VMobject()
            res_list = [Circle_Inversion(Dot(self.input_to_graph_point(x, grap)), about_point, inversion_radius).get_center() for x in np.arange(x_min, x_max, step_size)]
            res_vm.set_points_as_corners(res_list)
            return res_vm
        rad_of_inv = 2
        Inv_Circle = Circle(radius = rad_of_inv).set_color(PURPLE)
        Orig = Dot().scale(0.75)
        Inv_center = TexMobject("O").scale(0.75)
        Inv_center.next_to(Orig, direction = DL, buff = 0.00)
        parabola = self.get_graph(
            lambda x: x ** 2 + 1e-8,
            x_min = -10,
            x_max = 10,
            color = BLUE,
        )
        parabola_inv = Circle_Inverse_Func(parabola, inversion_radius = rad_of_inv)
        self.play(Write(Inv_Circle))
        self.wait(2)
        self.play(FadeInFromLarge(parabola))
        self.wait(2)
        parab_d = Dot(self.input_to_graph_point(-10, parabola)).scale(0.75)
        self.play(
            MoveAlongPath(parab_d, parabola),
            ShowCreation(parabola_inv),
            run_time = 10,
            rate_func = linear
        )
        self.wait(2)
        self.play(FadeOut(parabola), FadeOut(parabola_inv), FadeOut(parab_d))
        self.wait(2)
        self.camera_frame.save_state()
        title_card = TexMobject("\\underline{\\text{Generalized circles invert to Generalized circles}}")
        title_card.move_to(4 * UP)
        self.play(Write(title_card), self.camera_frame.shift, 0.5 * UP)
        self.wait(2)
        c = Circle(arc_center = 0.5 * RIGHT, radius = 0.5)
        self.play(Write(Orig), Write(Inv_center), Write(c))
        self.wait(2)
        ptA = Dot(RIGHT).scale(0.75)
        labelA = TexMobject("A").scale(0.75).next_to(ptA, direction = DR, buff = 0)
        ptAd = Circle_Inversion(ptA, inversion_radius = rad_of_inv).scale(0.75)
        labelAd = TexMobject("A'").scale(0.75).next_to(ptAd, direction = DR, buff = 0)
        rayOA = Line(ORIGIN, ptAd.get_center()).set_color(GRAY)
        self.play(Write(ptA), Write(labelA))
        self.wait(2)
        self.play(Write(ptAd), Write(labelAd))
        self.wait(2)
        self.play(Write(rayOA))
        self.wait(2)
        inv = Circle_Inversion(c, inversion_radius = rad_of_inv)
        alp = 75 * DEGREES
        ptB = ptA.copy()
        ptB.rotate(alp, about_point = (ptA.get_center() + Orig.get_center()) / 2)
        labelB = TexMobject("B").scale(0.75).next_to(ptB, direction = UL, buff = 0)
        ptBd = Circle_Inversion(ptB, inversion_radius = rad_of_inv).scale(0.75)
        labelBd = TexMobject("B'").scale(0.75).next_to(ptBd, direction = DR, buff = 0)
        rayOB = Line(ORIGIN, ptBd.get_center()).set_color(GRAY)
        rayBA = DashedLine(ptB.get_center(), ptA.get_center()).set_color(GRAY)
        rayBdAd = DashedLine(ptBd.get_center(), ptAd.get_center()).set_color(GRAY)
        self.play(Write(ptB), Write(labelB))
        self.wait(2)
        self.play(Write(ptBd), Write(labelBd), Write(rayOB))
        self.wait(2)
        self.play(Write(rayBA), Write(rayBdAd))
        self.wait(2)
        c.rotate(alp)
        Bright = Elbow(angle = PI + rayOB.get_angle()).shift(ptB.get_center()).set_color(YELLOW)
        Adright = Elbow(angle = PI / 2).shift(ptAd.get_center()).set_color(YELLOW)
        self.play(Write(Bright))
        self.wait(2)
        text_B = VGroup(TexMobject("B\\text{ is an angle in a semicircle}"), TexMobject("\\therefore \\angle B \\text{ is right angled}"))
        text_B.scale(0.75)
        text_B.arrange(direction = DOWN, buff = 0.5)
        text_B[1].align_to(text_B[0], direction = LEFT)
        text_B.move_to(5 * LEFT + 2 * UP)
        self.play(Write(text_B), self.camera_frame.shift, LEFT, title_card.shift, LEFT)
        self.wait(2)
        text_tri_sim = VGroup(TexMobject("\\text{We know }\\triangle{OAB}\\equiv\\triangle{OB'A'}"), TexMobject("\\implies \\angle OBA = \\angle OA'B'"), TexMobject("\\therefore \\angle OA'B' \\text{ is right angled}"))
        text_tri_sim.scale(0.75)
        text_tri_sim.arrange(direction = DOWN, buff = 0.5)
        text_tri_sim.move_to(1 * DOWN)
        for obj in text_tri_sim:
            obj.align_to(text_B[0], direction = LEFT)
        self.play(Write(text_tri_sim))
        self.wait(2)
        self.play(Write(Adright))
        self.wait(2)
        B_vg = VGroup(rayOB, rayBA, rayBdAd, Bright, ptBd, ptB, labelB, labelBd)
        def pt_updater(obj):
            br, bar, bdadr, bright, bd, b, lb, lbd = obj
            lb.next_to(b, direction = UL, buff = 0)
            new_bd = Circle_Inversion(b, inversion_radius = rad_of_inv).scale(0.75)
            lbd.next_to(new_bd, direction = RIGHT, buff = 0)
            bd.become(new_bd)
            br.become(Line(ORIGIN, bd.get_center()).set_color(GRAY))
            bar.become(DashedLine(b.get_center(), ptA.get_center()).set_color(GRAY))
            bdadr.become(DashedLine(bd.get_center(), ptAd.get_center()).set_color(GRAY))
            bright.become(Elbow(angle = PI + br.get_angle()).shift(b.get_center()).set_color(YELLOW))
        B_vg.add_updater(pt_updater)
        self.add(B_vg)
        self.play(MoveAlongPath(ptB, Arc(arc_center = 0.5 * RIGHT, radius = 0.5, start_angle = alp, angle = -alp / 2)), run_time = 5, rate_func = there_and_back)
        B_vg.clear_updaters()
        self.wait(2)
        self.play(Write(inv))
        self.wait(2)
        self.play(*[FadeOut(i) for i in [text_B, text_tri_sim, Bright, Adright, inv, rayBA, rayBdAd, rayOA, rayOB]], self.camera_frame.shift, 3 * RIGHT, title_card.shift, 3 * RIGHT)
        self.wait(2)
        c.rotate(-alp)
        A_grp = VGroup(ptA, ptAd, labelA, labelAd)
        B_grp = VGroup(ptB, ptBd, labelB, labelBd)
        c_grp = VGroup(c, A_grp, B_grp)
        def cir_updater(obj):
            cir, ag, bg = obj
            ag[0].become(Dot(cir.get_start()).scale(0.75))
            new_pb = Dot(cir.get_start()).scale(0.75)
            new_pb.rotate(alp, about_point = cir.get_center())
            bg[0].become(new_pb)
            def point_updater(obj, d = DR, dd = DR):
                a, ad, la, lad = obj
                la.next_to(a, direction = d, buff = 0)
                new_ad = Circle_Inversion(a, inversion_radius = rad_of_inv).scale(0.75)
                ad.become(new_ad)
                lad.next_to(ad, direction = dd, buff = 0)
                return VGroup(a, ad, la, lad)
            new_ag = point_updater(ag)
            new_bg = point_updater(bg, d = UL, dd = UR)
            ag.become(new_ag)
            bg.become(new_bg)
        c_grp.add_updater(cir_updater)
        self.add(c_grp)
        self.play(c.shift, 0.75 * RIGHT)
        self.wait(2)
        c_grp.clear_updaters()
        self.play(labelA.next_to, ptA, {"direction": DL, "buff": 0})
        self.wait(2)
        ptC = Dot(0.75 * RIGHT).scale(0.75)
        labelC = TexMobject("C").scale(0.75).next_to(ptC, direction = DL, buff = 0)
        ptCd = Circle_Inversion(ptC, inversion_radius = rad_of_inv).scale(0.75)
        labelCd = TexMobject("C'").scale(0.75).next_to(ptCd, direction = DL, buff = 0)
        rayOC = Line(ORIGIN, ptCd.get_center()).set_color(GRAY)
        rayBA = DashedLine(ptB.get_center(), ptA.get_center()).set_color(GRAY)
        rayBdAd = DashedLine(ptAd.get_center(), ptBd.get_center()).set_color(GRAY)
        rayOB = Line(ORIGIN, ptBd.get_center()).set_color(GRAY)
        rayBC = DashedLine(ptB.get_center(), ptC.get_center()).set_color(GRAY)
        rayBdCd = DashedLine(ptBd.get_center(), ptCd.get_center()).set_color(GRAY)
        C_grp = VGroup(ptC, ptCd, labelC, labelCd)
        self.play(Write(C_grp))
        self.wait(2)
        self.play(*[Write(i) for i in [rayOC, rayBA, rayBdAd, rayOB, rayBC, rayBdCd]])
        self.wait(2)
        self.play(self.camera_frame.scale, 0.5, self.camera_frame.move_to, ptCd.get_center() / 2)
        angleOBA = Sector(
            arc_center = ptB.get_center(),
            outer_radius = 0.25,
            start_angle = PI + rayOB.get_angle(),
            angle = PI + rayBA.get_angle() - rayOB.get_angle(),
            fill_color = YELLOW)
        self.play(Write(angleOBA))
        text_BeqAd = TexMobject("\\because \\triangle OBA \\equiv \\triangle OA'B'").scale(0.3125).move_to(4 * RIGHT + 1 * DOWN)
        angleOAdBd = Sector(
            arc_center = ptAd.get_center(),
            outer_radius = 0.25,
            start_angle = rayBdAd.get_angle(),
            angle = PI - rayBdAd.get_angle(),
            fill_color = YELLOW)
        self.play(Write(angleOAdBd))
        self.wait(2)
        angleABBd = Sector(
            arc_center = ptB.get_center(),
            outer_radius = 0.25,
            start_angle = rayBA.get_angle(),
            angle = - rayBA.get_angle() + rayOB.get_angle(),
            fill_color = GREEN)
        angleBdAdCd = Sector(
            arc_center = ptAd.get_center(),
            outer_radius = 0.25,
            start_angle = 0,
            angle = rayBdAd.get_angle(),
            fill_color = GREEN)
        self.play(Write(angleABBd), Write(angleBdAdCd))
        self.wait(2)
        self.play(FadeOut(angleOBA), FadeOut(angleOAdBd))
        self.wait(2)
        angleOBC = Sector(
            arc_center = ptB.get_center(),
            outer_radius = 0.25,
            start_angle = PI + rayOB.get_angle(),
            angle = PI + rayBC.get_angle() - rayOB.get_angle(),
            fill_color = BLUE)
        self.play(Write(angleOBC))
        self.wait(2)
        angleOCdBd = Sector(
            arc_center = ptCd.get_center(),
            outer_radius = 0.25,
            start_angle = PI,
            angle = rayBdCd.get_angle(),
            fill_color = BLUE)
        self.play(Write(angleOCdBd))
        self.wait(2)
        angleCBA = Elbow(angle = rayBA.get_angle() - PI / 2).shift(ptB.get_center()).set_color(YELLOW)
        angleCdBdAd = Elbow(angle = rayBdAd.get_angle() + PI).shift(ptBd.get_center()).set_color(YELLOW)
        self.play(Write(angleCBA), Write(angleCdBdAd))
        self.wait(2)
        self.play(
            *[FadeOut(i) for i in [angleOBC, angleABBd, angleBdAdCd, angleOCdBd]],
            self.camera_frame.scale, 2,
            self.camera_frame.move_to, 2 * RIGHT,
            title_card.move_to, 3.5 * UP + 2 * RIGHT)
        self.wait(2)
        bang = Line(c.get_center(), ptB.get_center()).get_angle()
        bdang = Line((ptAd.get_center() + ptCd.get_center()) / 2, ptBd.get_center()).get_angle()
        barc = Arc(arc_center = c.get_center(), radius = 0.5, start_angle = bang, angle = PI - bang)
        inv_barc = Circle_Inversion(barc, inversion_radius = rad_of_inv)
        bcirc_grp = VGroup(rayOB, rayBC, rayBA, rayBdAd, rayBdCd, ptBd, ptB, labelB, labelBd, angleCBA, angleCdBdAd)
        def invcircgrp(obj):
            br, bc, ba, bdad, bdcd, bd, b, lb, lbd, cba, cdbdad = obj
            new_bd = Circle_Inversion(b, inversion_radius = rad_of_inv).scale(0.75)
            bd.become(new_bd)
            lb.next_to(b, direction = UL, buff = 0)
            lbd.next_to(bd, direction = UR, buff = 0)
            br.become(Line(ORIGIN, bd.get_center()).set_color(GRAY))
            bc.become(DashedLine(b.get_center(), ptC.get_center()).set_color(GRAY))
            ba.become(DashedLine(b.get_center(), ptA.get_center()).set_color(GRAY))
            bdad.become(DashedLine(bd.get_center(), ptAd.get_center()).set_color(GRAY))
            bdcd.become(DashedLine(bd.get_center(), ptCd.get_center()).set_color(GRAY))
            cba.become(Elbow(angle = ba.get_angle() - PI / 2).shift(b.get_center()).set_color(YELLOW))
            cdbdad.become(Elbow(angle = bdad.get_angle()).shift(bd.get_center()).set_color(YELLOW))
        bcirc_grp.add_updater(invcircgrp)
        self.add(bcirc_grp)
        self.play(MoveAlongPath(ptB, barc), ShowCreation(inv_barc), run_time = 3, rate_func = linear)
        self.play(MoveAlongPath(ptB, barc), run_time = 3, rate_func = lambda t: linear(1 - t))
        bcirc_grp.clear_updaters()
        self.wait(2)

class InversionPreserves(GraphScene, MovingCameraScene):
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
        "default_graph_colors": [YELLOW, GREEN, YELLOW],
        "default_derivative_color": GREEN,
        "default_input_color": YELLOW,
        "default_riemann_start_color": YELLOW,
        "default_riemann_end_color": GREEN,
        "area_opacity": 0.8,
        "num_rects": 50,
        "dot_kwargs": {
            "radius": 0.25,
        },
        "line_kwargs": {
            "stroke_width": 2,
        },
        "fill_triangle_kwargs": {
            # "fill_color": YELLOW,
            "fill_opacity": .5,
            "stroke_width": 0,
        },
    }

    def setup(self):
        GraphScene.setup(self)
        MovingCameraScene.setup(self)
    
    def construct(self):
        self.setup_axes(animate = False)
        def Circle_Inverse_Func(grap, x_min = -5, x_max = 5, step_size = 0.01, about_point = ORIGIN, inversion_radius = 1):
            res_vm = VMobject()
            res_list = [Circle_Inversion(Dot(self.input_to_graph_point(x, grap)), about_point, inversion_radius).get_center() for x in np.arange(x_min, x_max, step_size)]
            res_vm.set_points_as_corners(res_list)
            return res_vm
        rad_of_inv = 2
        Inv_Circle = Circle(radius = rad_of_inv).set_color(PURPLE)
        Orig = Dot().scale(0.75)
        Inv_center = TexMobject("O").scale(0.75)
        Inv_center.next_to(Orig, direction = DL, buff = 0.00)
        parabola = self.get_graph(
            lambda x: x ** 3 + 1e-8,
            x_min = -5,
            x_max = 5,
            color = BLUE,
        )
        parabola_inv = Circle_Inverse_Func(parabola, inversion_radius = rad_of_inv)
        self.play(Write(Inv_Circle))
        self.wait(2)
        self.play(Write(parabola), Write(parabola_inv))
        self.wait(2)
        self.play(FadeOut(parabola), FadeOut(parabola_inv))
        self.wait(2)
        title_card = TexMobject("\\underline{\\text{Invariants under Inversion}}")
        title_card.move_to(3.5 * UP)
        self.play(Write(title_card))
        ptP = Dot().scale(0.75).move_to(5 * RIGHT)
        labelP = TexMobject("P").scale(0.75).next_to(ptP, direction = DR, buff = 0)
        t_lines = CircleTangentLines(Inv_Circle, ptP)
        self.play(Write(ptP), Write(labelP))
        self.wait(2)
        t_points = VGroup()
        for t in t_lines:
            t_points.add(Dot().scale(0.75).move_to(t.get_end()))
        self.play(Write(t_lines), Write(t_points))
        self.wait(2)
        ortho_circ = Circle(arc_center = ptP.get_center(), radius = t_lines[0].get_length())
        self.play(Write(ortho_circ))
        self.wait(2)
        scaled_t_lines = VGroup()
        for t in t_points:
            scaled_t_lines.add(DashedLine(ptP.get_center(), ptP.get_center() + 1.25 * (t.get_center() - ptP.get_center())).set_color(GREY))
        radial_lines = VGroup()
        for t in t_points:
            radial_lines.add(DashedLine(ORIGIN, 1.25 * t.get_center()).set_color(GREY))
        self.play(Write(scaled_t_lines), Write(radial_lines))
        self.wait(2)
        r_angles = VGroup()
        r_angles.add(Elbow(angle = PI / 2 + t_lines[0].get_angle()).shift(t_points[0].get_center()).set_color(YELLOW))
        r_angles.add(Elbow(angle = 2 * PI / 2 + t_lines[1].get_angle()).shift(t_points[1].get_center()).set_color(YELLOW))
        self.play(Write(r_angles))
        self.wait(2)
        self.play(*[FadeOut(i) for i in [labelP, t_lines, scaled_t_lines[1], t_points[1], r_angles[1], radial_lines[1]]])
        self.wait(2)
        ortho_grp = VGroup(radial_lines[0], scaled_t_lines[0], t_points[0], ortho_circ, r_angles[0], ptP)
        def ortho_updater(obj):
            rl, tl, tp, ort, ang, pp = obj
            new_tl = CircleTangentLines(Inv_Circle, pp)[0]
            tp.become(Dot().scale(0.75).move_to(new_tl.get_end()))
            tl.become(DashedLine(pp.get_center(), pp.get_center() + 1.25 * (tp.get_center() - pp.get_center())).set_color(GREY))
            rl.become(DashedLine(ORIGIN, 1.25 * tp.get_center()).set_color(GREY))
            ort.become(Circle(arc_center = pp.get_center(), radius = new_tl.get_length()))
            ang.become(Elbow(angle = PI / 2 + new_tl.get_angle()).shift(tp.get_center()).set_color(YELLOW))
        ortho_grp.add_updater(ortho_updater)
        self.add(ortho_grp)
        self.play(MoveAlongPath(ptP, Ellipse(width = 10, height = 5)), run_time = 5, rate_func = linear)
        self.wait(2)
        ortho_grp.clear_updaters()
        self.play(FadeOut(ortho_grp))
        self.wait(2)
        ell = Ellipse(width = 8, height = 4)
        inv_ell = Circle_Inversion(ell, inversion_radius = rad_of_inv)
        self.play(Write(ell), Write(inv_ell))
        temp_circ = Circle().move_to(6 * RIGHT).set_color(GREEN)
        inv_temp_circ = Circle_Inversion(temp_circ, inversion_radius = rad_of_inv)
        temp_c_grp = VGroup(temp_circ, inv_temp_circ)
        def inv_upd(obj):
            c, ic = obj
            ic.become(Circle_Inversion(c, inversion_radius = rad_of_inv))
        temp_c_grp.add_updater(inv_upd)
        self.add(temp_c_grp)
        self.play(temp_circ.shift, LEFT, run_time = 2)
        self.wait(2)
        self.play(temp_circ.shift, LEFT, run_time = 2)
        self.wait(2)
        self.play(temp_circ.shift, 1.25 * LEFT, run_time = 2)
        self.wait(2)
        self.play(temp_circ.move_to, 6 * LEFT, run_time = 2)
        self.wait(2)
        temp_c_grp.clear_updaters()
        self.wait(2)
        self.play(*[FadeOut(i) for i in [temp_c_grp, ell, inv_ell]])
        self.wait(2)
        circ1 = Circle().move_to(RIGHT)
        inv_circ1 = Circle_Inversion(circ1, inversion_radius = rad_of_inv)
        circ2 = Circle().move_to(UP)
        inv_circ2 = Circle_Inversion(circ2, inversion_radius = rad_of_inv)
        self.play(Transform(title_card, TexMobject("\\underline{\\text{Where is the Second intersection point?}}").move_to(3.5 * UP)))
        self.play(Write(Orig), Write(Inv_center))
        self.wait(2)
        self.play(Write(circ1))
        self.wait(2)
        self.play(Write(inv_circ1))
        self.wait(2)
        self.play(Write(circ2), Write(inv_circ2))
        self.wait(2)

if __name__ == "__main__":
    module_name = os.path.basename(__file__)
    #command_A = "manim -p  --video_dir ~/Downloads/  "
    clear_cmd = "cls"
    command_A = "python -m manim "
    #command_B = module_name + " " + "IntroductionScene" + " -pl -n 51,55"
    #command_B = module_name + " " + "IntroductionScene" + " -p"
    command_B = module_name + " -p"
    os.system(clear_cmd)
    os.system(command_A + command_B)