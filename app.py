import numpy as np
import pandas as pd
import plotly.graph_objects as go
import streamlit as st

st.set_page_config(
    page_title="Навчальний симулятор rAAV", page_icon="🧬", layout="wide"
)

st.title("🧬 Інтерактивний симулятор біопроцесу rAAV (HEK293)")
st.caption(
    "Модель upstream-виробництва рекомбінантного аденоасоційованого вірусу."
)

# --- ТАБИ НАВІГАЦІЇ ---
tab1, tab2, tab3, tab4 = st.tabs([
    "🎮 Симулятор та розрахунки",
    "🖼️ Візуалізація біореактора",
    "📚 Картки-шпори (Покращення якості)",
    "⚖️ Масштабування (200L ➔ 2000L)",
])

# --- SIDEBAR: ВХІДНІ ПАРАМЕТРИ ---
st.sidebar.header("⚙️ Параметри біопроцесу")

cell_density = st.sidebar.slider(
    "Живі клітини при трансфекції (Xₜ), 10⁶ cells/mL",
    0.5,
    8.0,
    3.0,
    0.1,
    help="Оптимум: 2.5–3.5 x 10⁶ cells/mL",
)

dna_per_cell = st.sidebar.slider(
    "Загальна ДНК на клітину (Dₜ), ум. од.",
    0.2,
    3.0,
    1.0,
    0.1,
    help="Передозування ДНК призводить до появи порожніх капсидів",
)

st.sidebar.markdown("**Співвідношення плазмід (Сума = 1.0)**")
ratio_transgene = st.sidebar.slider(
    "Касета трансгена (GoI)", 0.1, 0.7, 0.33, 0.01
)
ratio_repcap = st.sidebar.slider("Rep/Cap плазміда", 0.1, 0.7, 0.33, 0.01)
ratio_helper = st.sidebar.slider("Helper плазміда", 0.1, 0.7, 0.34, 0.01)

harvest_hpt = st.sidebar.slider(
    "Час збору (Harvest time), hpt",
    24,
    120,
    60,
    1,
    help="Години після трансфекції (hours post-transfection)",
)

serotype = st.sidebar.selectbox(
    "Серотип rAAV", options=["AAV2", "AAV5", "AAV8", "AAV9"]
)

lysis_mode = st.sidebar.selectbox(
    "Стратегія лізису",
    options=[
        "Без лізису (збір супернатанту)",
        "Частковий лізис",
        "Повний лізис",
    ],
)

scale_l = st.sidebar.select_slider(
    "Масштаб біореактора, L", options=[2, 50, 200, 500, 1000, 2000], value=200
)

# --- РОЗРАХУНКОВИЙ БЛОК ---
base_titer = 3.0e10
cell_effect = max(0.1, 1.0 - abs(cell_density - 3.0) / 8.0)
dna_effect = max(
    0.1, 1.0 + 0.25 * (dna_per_cell - 1.0) - 0.20 * max(0.0, dna_per_cell - 1.5)
)

plasmid_sum = ratio_transgene + ratio_repcap + ratio_helper
imbalance = abs(plasmid_sum - 1.0) + abs(ratio_transgene - 0.33) * 0.5
balance_effect = max(0.1, 1.0 - 1.5 * imbalance)

if harvest_hpt < 48:
  time_effect = 0.6 + (harvest_hpt - 24) / 24 * 0.4
elif harvest_hpt <= 72:
  time_effect = 1.0
else:
  time_effect = 1.0 - (harvest_hpt - 72) / 48 * 0.25

intracellular_fraction = {
    "AAV2": 0.85,
    "AAV5": 0.60,
    "AAV8": 0.75,
    "AAV9": 0.70,
}[serotype]

lysis_release = {
    "Без лізису (збір супернатанту)": 0.0,
    "Частковий лізис": 0.55,
    "Повний лізис": 0.95,
}[lysis_mode]

available_fraction = (1.0 - intracellular_fraction) + (
    intracellular_fraction * lysis_release
)

scale_penalty = {2: 1.00, 50: 0.98, 200: 0.95, 500: 0.90, 1000: 0.85, 2000: 0.78}[
    scale_l
]

titer = (
    base_titer
    * cell_effect
    * dna_effect
    * balance_effect
    * time_effect
    * available_fraction
    * scale_penalty
)

base_full = 0.28
full_fraction = (
    base_full
    - 0.08 * max(0.0, dna_per_cell - 1.0)
    - 0.15 * imbalance
    - 0.05 * max(0.0, harvest_hpt - 72) / 24.0
    - 0.03 * max(0.0, cell_density - 5.0)
)
full_pct = max(5.0, min(45.0, full_fraction * 100))

impurity_load = (
    1.0
    + 2.0 * lysis_release
    + 0.5 * max(0.0, harvest_hpt - 72) / 24.0
    + 0.3 * max(0.0, cell_density - 4.0)
)

quality_score = (
    0.50 * (full_pct / 30.0 * 100)
    + 0.30 * max(0.0, (1.0 - (impurity_load - 1.0) / 3.0) * 100)
    + 0.20 * min(100.0, (titer / 3.0e10) * 100)
)
quality_score = max(0.0, min(100.0, quality_score))

harvest_volume_ml = scale_l * 1000.0 * 0.9
total_vg = titer * harvest_volume_ml

# --- ТАБ 1: СИМУЛЯТОР ТА РОЗРАХУНКИ ---
with tab1:
  st.subheader("📊 Головні показники лоту")
  c1, c2, c3, c4 = st.columns(4)
  c1.metric("Геномний титр (Crude)", f"{titer:.2e} vg/mL")
  c2.metric("Повний вихід", f"{total_vg:.2e} vg")
  c3.metric("Частка повних капсидів", f"{full_pct:.1f} %")
  c4.metric("Індекс якості лоту", f"{quality_score:.0f} / 100")

  st.markdown("---")
  st.subheader("⚠️ Попередження та критичні компроміси")

  if abs(plasmid_sum - 1.0) > 0.02:
    st.error(f"❌ Сума часток плазмід має дорівнювати 1.00 (Зараз: {plasmid_sum:.2f})")

  if titer > 3.0e10 and full_pct < 20.0:
    st.warning(
        "⚡ **Високий титр, але низька якість!** Багато часток є порожніми"
        " (без генома). Геномний титр ≠ інфекційний титр."
    )

  if dna_per_cell > 1.8:
    st.warning(
        "🧬 **Надлишок ДНК на клітину:** знижує частку повних капсидів та"
        " підвищує цитотоксичність."
    )

  if harvest_hpt > 84:
    st.warning(
        "⏱️ **Пізній збір (>84 hpt):** деградація вірусу та масовий накопичення"
        " ДНК/білків хазяїна (HCDNA/HCP)."
    )

  if lysis_mode == "Повний лізис":
    st.warning(
        "🧪 **Повний лізис:** збільшує титр, але критично забруднює лізат"
        " клітинним сміттям (складніше downstream-очищення)."
    )

  if scale_l >= 1000:
    st.info(
        "🏭 **Промисловий масштаб (1000-2000 L):** фізико-хімічні обмеження"
        " (mixing time, kLa, shear stress) знижують питому продуктивність."
    )

  # Динамічний графік кінетики
  st.subheader("📈 Динаміка біопроцесу (0 ➔ Harvest)")
  t_steps = np.linspace(0, harvest_hpt, 100)
  xv_t = cell_density * np.exp(0.015 * t_steps) * (1 - t_steps / 150)
  vg_t = (
      (titer / (1 + np.exp(-0.1 * (t_steps - 36))))
      if harvest_hpt > 24
      else np.zeros(100)
  )

  fig = go.Figure()
  fig.add_trace(
      go.Scatter(
          x=t_steps,
          y=xv_t,
          name="Щільність клітин (10⁶/mL)",
          line=dict(color="blue", width=3),
      )
  )
  fig.add_trace(
      go.Scatter(
          x=t_steps,
          y=vg_t,
          name="Титр AAV (vg/mL)",
          yaxis="y2",
          line=dict(color="green", width=3, dash="dash"),
      )
  )

  fig.update_layout(
      title="Кінетика росту клітин та накопичення rAAV у біореакторі",
      xaxis=dict(title="Час після трансфекції (hpt), годин"),
      yaxis=dict(title="Живі клітини (10⁶ cells/mL)", side="left"),
      yaxis2=dict(
          title="Титр rAAV (vg/mL)", overlaying="y", side="right", type="log"
      ),
      legend=dict(orientation="h", y=1.15),
      height=400,
  )
  st.plotly_chart(fig, use_container_width=True)

# --- ТАБ 2: ВІЗУАЛІЗАЦІЯ БІОРЕАКТОРА ---
with tab2:
  st.subheader("🖼️ Інтерактивна схема Stirred-Tank Біореактора")

  col_vis, col_status = st.columns([1, 1])

  with col_vis:
    # SVG схема біореактора
    agitator_speed = "Швидка" if scale_l >= 1000 else "Оптимальна"
    shear_warning = scale_l >= 1000
    reactor_color = "#e3f2fd" if not shear_warning else "#fff3e0"

    svg_code = f"""
        <svg width="350" height="420" viewBox="0 0 350 420" xmlns="http://www.w3.org/2000/svg">
          <!-- Bioreactor Tank Body -->
          <rect x="75" y="60" width="200" height="280" rx="30" ry="30" fill="{reactor_color}" stroke="#1565c0" stroke-width="6"/>
          <!-- Fluid Level -->
          <rect x="80" y="120" width="190" height="215" rx="10" ry="10" fill="#90caf9" opacity="0.6"/>
          <!-- Agitator Shaft -->
          <line x1="175" y1="30" x2="175" y2="300" stroke="#424242" stroke-width="8"/>
          <!-- Impeller Blades -->
          <rect x="125" y="240" width="100" height="15" fill="#616161" rx="4"/>
          <rect x="135" y="280" width="80" height="12" fill="#616161" rx="4"/>
          <!-- Motor -->
          <rect x="145" y="10" width="60" height="30" fill="#212121" rx="5"/>
          <!-- Bubbles / Sparger -->
          <circle cx="120" cy="300" r="4" fill="#ffffff" opacity="0.8"/>
          <circle cx="150" cy="270" r="6" fill="#ffffff" opacity="0.8"/>
          <circle cx="190" cy="290" r="5" fill="#ffffff" opacity="0.8"/>
          <circle cx="220" cy="260" r="4" fill="#ffffff" opacity="0.8"/>
          <!-- Labels -->
          <text x="175" y="100" font-family="Arial" font-size="14" font-weight="bold" text-anchor="middle" fill="#0d47a1">Об'єм: {scale_l} L</text>
          <text x="175" y="170" font-family="Arial" font-size="12" text-anchor="middle" fill="#1b5e20">HEK293: {cell_density}M/mL</text>
          <text x="175" y="190" font-family="Arial" font-size="12" text-anchor="middle" fill="#b71c1c">Серотип: {serotype}</text>
        </svg>
        """
    st.components.v1.html(svg_code, height=430)

  with col_status:
    st.markdown("### 📋 Фізико-хімічний стан біосистеми")
    st.write(f"**Об'єм реактора:** `{scale_l} L`")
    st.write(f"**Швидкість змішування:** `{agitator_speed}`")
    st.write(
        "**Температура:** `37.0 °C` | **pH:** `7.1` | **Розчинений кисень"
        " (DO):** `40%`"
    )

    if scale_l >= 1000:
      st.error(
          "🔴 **Високе зсувне напруження (Shear Stress):** Клітини HEK293"
          " зазнають механічного пошкодження від лопатей мішалки."
      )
      st.warning(
          "🟡 **Обмеження kLa (Масоперенос кисню):** Потрібне додаткове"
          " збагачення чистого кисню."
      )
    else:
      st.success(
          "🟢 **Оптимальний гомеостаз:** Змішування та киснепостачання"
          " рівномірні по всьому об'єму."
      )

    st.markdown("---")
    st.markdown("**Розподіл продукту (Intracellular vs Extracellular):**")
    st.progress(intracellular_fraction)
    st.caption(
        f"Внутрішньоклітинний AAV: {intracellular_fraction*100:.0f}% | У"
        f" супернатанті: {(1-intracellular_fraction)*100:.0f}%"
    )

# --- ТАБ 3: КАРТКИ-ШПОРИ ---
with tab3:
  st.subheader("📚 Картки-шпори біоінженера: Як досягти максимальної якості")

  col_card1, col_card2 = st.columns(2)

  with col_card1:
    st.info("""
        ### 🎯 1. Як збільшити частку повних капсидів (% Full)
        
        * **Баланс плазмід:** Тримайте молярне співвідношення `Transgene : Rep/Cap : Helper` близьким до **1:1:1**. Дисбаланс призводить до синтезу "порожніх оболонок".
        * **Контроль дози ДНК:** Не перевищуйте **1.0–1.2 мкг ДНК / 10⁶ клітин**. Надлишок ДНК викличує цитотоксичність та утворення неінфекційних часткових часток.
        * **Оптимальне вікно збору:** Збирайте біомасу строго у проміжку **48–72 hpt**. Після 72 hpt частка повних капсидів починає швидко падати через деградацію.
        """)

    st.success("""
        ### 🛡️ 3. Як оптимізувати трансфекцію HEK293
        
        * **Щільність при трансфекції (Xₜ):** Вводьте комплекс строго при **2.5–3.5 × 10⁶ cells/mL** у фазі експоненційного росту.
        * **Час формування комплексу DNA:PEI:** Комплекс має формуватися **10–15 хвилин**. Подача понад 20 хвилин знижує ефективність трансфекції у 2–3 рази!
        """)

  with col_card2:
    st.warning("""
        ### 🧪 2. Як мінімізувати домішки (HCDNA / HCP)
        
        * **Відмова від агресивного лізису:** Використовуйте **Частковий лізис** або збір супернатанту (якщо серотип це дозволяє, наприклад AAV5), щоб не руйнувати ядра клітин.
        * **Ферментативна обробка (Benzonase):** Додавайте ендонуклеазу одразу після збору для розщеплення в'язкої хазяйської ДНК.
        * **Зсувний контроль (Shear control):** Знижуйте оберти імпелера перед збором, щоб запобігти передчасному механічному некрозу клітин.
        """)

    st.error("""
        ### 🏭 4. Чек-лист успішного Scale-Up (2000 L)
        
        * **P/V & kLa:** Масштабуйте процес за рівною питомою потужністю перемішування ($P/V$) та коефіцієнтом масопереносу кисню ($k_L a$).
        * **Inline-трансфекція:** Застосовуйте статичні міксери для швидкої автоматичної подачі комплексу ДНК у великий біореактор.
        * **Downstream готовий:** Будьте готові до втрат 30–50% продукту на стадії IEX-хроматографії при розділенні повних/порожніх капсидів.
        """)

# --- ТАБ 4: МАСШТАБУВАННЯ ---
with tab4:
  st.subheader("⚖️ Порівняльний аналіз масштабування: 200 L ➔ 2000 L")

  case_data = {
      "Показник": [
          "Живі клітини при трансфекції",
          "Час збору",
          "Прогнозований титр (vg/mL)",
          "Частка повних капсидів (% full)",
          "Навантаження домішками",
          "Критичне лімітування",
      ],
      "Пілотний біореактор (200 L)": [
          "3.0 x 10⁶ cells/mL",
          f"{harvest_hpt} hpt",
          f"{titer:.2e}",
          f"{full_pct:.1f}%",
          "Низьке / Базове",
          "Швидке однорідне змішування",
      ],
      "Промисловий біореактор (2000 L)": [
          "3.0 x 10⁶ cells/mL",
          f"{harvest_hpt} hpt",
          f"{titer * (0.78 / scale_penalty):.2e}",
          f"{max(5.0, full_pct - 3.0):.1f}%",
          "Підвищене",
          "Час подачі комплексу, Shear stress, kLa",
      ],
  }
  st.dataframe(pd.DataFrame(case_data), use_container_width=True)

  st.success(
      "💡 **Ключове інженерне резюме:** Головний критерій успішного scale-up —"
      " збереження питомого виходу повних функціональних векторів на 1 літр"
      " середовища та збереження CQA (Critical Quality Attributes), а не"
      " максимальний титр будь-якою ціною."
  )
